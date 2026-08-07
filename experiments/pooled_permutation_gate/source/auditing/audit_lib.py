#!/usr/bin/env python3
"""Attribution and auditing library for a fitted Decision Mesh.

Everything here is exact linear algebra CONDITIONAL on three selected
objects: the mesh, the converged working weights, and the shrinkage
factors. See README.md for what that conditioning means and for the two
stated approximations (reconstructed weights; staged-IRLS vs one-shot
WLS heights).

State expected: <prefix>_mesh.csv, <prefix>_hier_vertices.csv, the pool
csv (wala, wac, n0, k_term), <prefix>_obs.csv, and law_shared.json.
"""
import json
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import matplotlib.tri as mtri

STRATA_EDGES = [128, 512, 4096]


class MeshAudit:
    def __init__(self, pools_csv, prefix, law_json):
        d = pd.read_csv(pools_csv)
        self.n = d.n0.values.astype(float)
        self.k = d.k_term.values.astype(float)
        self.wala, self.wac = d.wala.values, d.wac.values
        self.held = np.arange(len(d)) % 2 == 1
        self.eta0 = float(np.log(self.k.sum() / self.n.sum()
                                 / (1 - self.k.sum() / self.n.sum())))
        obs = pd.read_csv(f"{prefix}_obs.csv")
        ph = np.clip(obs.p_hat.values, 1e-9, 1 - 1e-9)
        self.fit_lo = np.log(ph / (1 - ph))

        # COLUMN SEMANTICS (pinned from mesh.cpp/vertex.cpp, July 22):
        # 'height' is the COMMITTED production value (equals the mesh dump and
        # p_hat); 'new_height' is the last raw sweep proposal (pre-shrink).
        # An earlier version had these roles swapped; caught when the
        # 'shrunken' ramp profile dove below the data while p_hat did not.
        hv = pd.read_csv(f"{prefix}_hier_vertices.csv")
        hv = hv[hv.active == 1].reset_index(drop=True)
        self.hv = hv
        self.id2row = {int(i): r for r, i in enumerate(hv.id)}
        self.depth = hv.depth.values.astype(int)
        self.lam = hv.lambda_v.values
        self.admitted = hv.gate_admitted.values.astype(bool)
        # surpluses of the SHRUNKEN surface (what the model outputs)
        par = hv[["parent0", "parent1"]].values.astype(int)
        nh = hv.height.values  # production
        mid = np.where(par[:, 0] >= 0,
                       0.5 * (nh[[self.id2row.get(p, 0) for p in par[:, 0]]]
                              + nh[[self.id2row.get(p, 0) for p in par[:, 1]]]),
                       0.0)
        self.surplus = nh - np.where(par[:, 0] >= 0, mid, 0.0)
        self.par = par
        self.nh = nh
        # implied shrinkage factor per vertex: kept surplus / raw surplus.
        # (the dump's lambda_v column is a demoted diagnostic, not the
        # reshrink factor; verified by kernel localization failure)
        rh = hv.new_height.values  # raw sweep proposal
        mid_raw = np.where(par[:, 0] >= 0,
                           0.5 * (rh[[self.id2row.get(p, 0) for p in par[:, 0]]]
                                  + rh[[self.id2row.get(p, 0) for p in par[:, 1]]]),
                           0.0)
        raw_surplus = rh - np.where(par[:, 0] >= 0, mid_raw, 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            li = np.where(np.abs(raw_surplus) > 1e-9,
                          self.surplus / raw_surplus, 0.0)
        self.lam_implied = np.clip(li, 0.0, 1.0)

        # fine-mesh triangulation from the face dump (conforming; verified)
        mesh = pd.read_csv(f"{prefix}_mesh.csv")
        corners = np.vstack([mesh[[f"x{i}", f"y{i}", f"h{i}"]].values
                             for i in (0, 1, 2)])
        xy, inv = np.unique(np.round(corners[:, :2], 9), axis=0,
                            return_inverse=True)
        self.vxy = xy
        self.tris = inv.reshape(3, len(mesh)).T
        self.tri = mtri.Triangulation(xy[:, 0], xy[:, 1], triangles=self.tris)
        self.finder = self.tri.get_trifinder()
        # map triangulation vertices -> hier rows by coordinates
        key = {(round(x, 9), round(y, 9)): r
               for r, (x, y) in enumerate(zip(hv.x, hv.y))}
        self.tv2h = np.array([key[(round(x, 9), round(y, 9))]
                              for x, y in xy], dtype=int)
        self.face_area = 0.5 * np.abs(
            (mesh.x1 - mesh.x0) * (mesh.y2 - mesh.y0)
            - (mesh.x2 - mesh.x0) * (mesh.y1 - mesh.y0)).values

        # pool design: barycentric coords in the containing face
        self.ux = (self.wala - self.wala.min()) / (self.wala.max() - self.wala.min())
        self.uy = (self.wac - self.wac.min()) / (self.wac.max() - self.wac.min())
        self.pool_face = self.finder(np.clip(self.ux, 0, 0.999999),
                                     np.clip(self.uy, 0, 0.999999))
        law = json.load(open(law_json))["law"]
        tot = np.array([t[3] for t in law])
        sidx = np.searchsorted(STRATA_EDGES, self.n, side="right")
        self.s2 = tot[sidx]
        ph_e = np.clip((self.k + 0.5) / (self.n + 1), 1e-9, 1 - 1e-9)
        self.v_samp = 1.0 / ((self.n + 1) * ph_e * (1 - ph_e))
        # reconstructed working weights (logit-variance units; TRAIN half only,
        # mirroring the fit's held-out convention)
        self.w = np.where(self.held, 0.0, 1.0 / (self.v_samp + self.s2))
        self._build_design()

    # ---------------------------------------------------------------- design
    def _bary(self, fidx, x, y):
        t = self.tris[fidx]
        p = self.vxy[t]
        d = ((p[1, 1] - p[2, 1]) * (p[0, 0] - p[2, 0])
             + (p[2, 0] - p[1, 0]) * (p[0, 1] - p[2, 1]))
        a = ((p[1, 1] - p[2, 1]) * (x - p[2, 0])
             + (p[2, 0] - p[1, 0]) * (y - p[2, 1])) / d
        b = ((p[2, 1] - p[0, 1]) * (x - p[2, 0])
             + (p[0, 0] - p[2, 0]) * (y - p[2, 1])) / d
        return t, np.array([a, b, 1 - a - b])

    def _build_design(self):
        rows, cols, vals = [], [], []
        for i in np.where(self.pool_face >= 0)[0]:
            t, phi = self._bary(self.pool_face[i], self.ux[i], self.uy[i])
            rows += [i] * 3
            cols += list(t)
            vals += list(phi)
        nv = len(self.vxy)
        self.X = sp.csr_matrix((vals, (rows, cols)),
                               shape=(len(self.n), nv))
        W = sp.diags(self.w)
        self.M = (self.X.T @ W @ self.X).tocsc()
        self.M_solve = spl.factorized(self.M + 1e-10 * sp.identity(nv, format="csc"))

    # ---------------------------------------------------------------- ladder
    def ladder(self, wala_q, wac_q):
        """Exact hierarchical decomposition of the model value at a point.
        Returns rows (hier_row, depth, psi, surplus, contribution_logits)
        and the base rate. Sums exactly to the surface value."""
        xq = (wala_q - self.wala.min()) / (self.wala.max() - self.wala.min())
        yq = (wac_q - self.wac.min()) / (self.wac.max() - self.wac.min())
        f = int(self.finder(np.array([np.clip(xq, 0, 0.999999)]),
                            np.array([np.clip(yq, 0, 0.999999)]))[0])
        if f < 0:
            raise ValueError("query outside domain")
        t, phi = self._bary(f, xq, yq)
        acc = {}
        for tv, p in zip(t, phi):
            hr = self.tv2h[tv]
            acc[hr] = acc.get(hr, 0.0) + p
        out = []
        while acc:
            hr = max(acc, key=lambda r: self.depth[r])
            wgt = acc.pop(hr)
            out.append((hr, int(self.depth[hr]), wgt, self.surplus[hr],
                        wgt * self.surplus[hr] / 100.0))
            p0, p1 = self.par[hr]
            if p0 >= 0:
                for p in (p0, p1):
                    r = self.id2row[int(p)]
                    acc[r] = acc.get(r, 0.0) + wgt / 2.0
        return out, self.eta0

    def ladder_check(self, wala_q, wac_q):
        rows, eta0 = self.ladder(wala_q, wac_q)
        total = eta0 + sum(r[4] for r in rows)
        xq = (wala_q - self.wala.min()) / (self.wala.max() - self.wala.min())
        yq = (wac_q - self.wac.min()) / (self.wac.max() - self.wac.min())
        f = int(self.finder(np.array([xq]), np.array([yq]))[0])
        t, phi = self._bary(f, xq, yq)
        direct = eta0 + float(sum(phi[j] * self.nh[self.tv2h[t[j]]]
                                  for j in range(3))) / 100.0
        return total, direct

    # ---------------------------------------------------------------- kernel
    def _c_vector(self, wala_q, wac_q):
        """Coefficients of the model value on RAW (pre-shrink) fine-mesh
        heights, composing the shrinkage operator through the hierarchy."""
        xq = (wala_q - self.wala.min()) / (self.wala.max() - self.wala.min())
        yq = (wac_q - self.wac.min()) / (self.wac.max() - self.wac.min())
        f = int(self.finder(np.array([np.clip(xq, 0, 0.999999)]),
                            np.array([np.clip(yq, 0, 0.999999)]))[0])
        t, phi = self._bary(f, xq, yq)
        acc, craw = {}, {}
        for tv, p in zip(t, phi):
            acc[self.tv2h[tv]] = acc.get(self.tv2h[tv], 0.0) + p
        # Shrunken height: new_h = mid + lam*(raw_h - mid)
        #                        = lam*raw_h + (1-lam)*mid(parents).
        # So a weight wgt on a shrunken vertex places wgt*lam on its own raw
        # height and wgt*(1-lam)/2 on each parent's shrunken height (recurse).
        # Roots carry their raw height directly.
        while acc:
            hr = max(acc, key=lambda r: self.depth[r])
            wgt = acc.pop(hr)
            p0, p1 = self.par[hr]
            if p0 < 0:
                craw[hr] = craw.get(hr, 0.0) + wgt
            else:
                craw[hr] = craw.get(hr, 0.0) + wgt * self.lam_implied[hr]
                for p in (p0, p1):
                    r = self.id2row[int(p)]
                    acc[r] = acc.get(r, 0.0) + wgt * (1 - self.lam_implied[hr]) / 2.0
        return craw

    def kernel(self, wala_q, wac_q):
        """Per-pool influence weights of the estimate at (wala, wac).
        estimate ~= sum_i kernel_i * y_i with y in logits (+ eta0)."""
        craw = self._c_vector(wala_q, wac_q)
        c = np.zeros(len(self.vxy))
        h2tv = {h: tv for tv, h in enumerate(self.tv2h)}
        for hr, cv in craw.items():
            c[h2tv[hr]] = cv
        z = self.M_solve(c)
        kern = self.w * (self.X @ z)
        return kern

    def kernel_stats(self, kern):
        a = np.abs(kern)
        s = a.sum()
        # abs-based concentration (robust to sign cancellation)
        neff = (a.sum() ** 2) / np.maximum((kern ** 2).sum(), 1e-30)
        sidx = np.searchsorted(STRATA_EDGES, self.n, side="right")
        shares = [a[sidx == j].sum() / max(s, 1e-12) for j in range(4)]
        return {"n_eff": float(neff), "stratum_abs_shares": shares}

    def kernel_check(self, wala_q, wac_q):
        """eta0 + kernel . y  vs the actual surface value; the gap measures
        the stated staged-IRLS approximation."""
        kern = self.kernel(wala_q, wac_q)
        y = np.log(np.clip((self.k + 0.5) / (self.n + 1), 1e-9, 1 - 1e-9)
                   / np.clip(1 - (self.k + 0.5) / (self.n + 1), 1e-9, 1 - 1e-9)) - self.eta0
        pred = self.eta0 + float(kern @ y)
        _, direct = self.ladder_check(wala_q, wac_q)
        return pred, direct

    # ------------------------------------------------------------- influence
    def leverage(self):
        """Diagonal leverage per pool (train half), via dense inverse of the
        vertex normal matrix (fine at ~2k vertices)."""
        Minv = spl.inv(sp.csc_matrix(self.M + 1e-10 * sp.identity(
            self.M.shape[0], format="csc"))).toarray()
        lev = np.zeros(len(self.n))
        Xc = self.X.tocsr()
        for i in np.where(self.w > 0)[0]:
            r = Xc.getrow(i)
            idx, val = r.indices, r.data
            lev[i] = self.w[i] * val @ Minv[np.ix_(idx, idx)] @ val
        return lev

    def self_influence_map(self, lev):
        """|leave-one-out change of the surface at the pool's own location|,
        centi-logits: the mega-risk map."""
        resid = (np.log(np.clip((self.k + 0.5) / (self.n + 1), 1e-9, 1 - 1e-9)
                        / np.clip(1 - (self.k + 0.5) / (self.n + 1), 1e-9, 1 - 1e-9))
                 - self.fit_lo)
        # leverage ~ 1 means a sole-support pool: its removal leaves the
        # local surface unconstrained and the LOO delta is unbounded. Cap
        # the ratio and let callers report the sole-support count.
        capped = np.minimum(lev, 0.99)
        out = np.where(self.w > 0,
                       capped * np.abs(resid) / (1 - capped) * 100, 0.0)
        self.n_sole_support = int(((lev > 0.99) & (self.w > 0)).sum())
        return out

    # -------------------------------------------------------------- ablation
    def psi_matrix(self, rows_subset):
        """Hierarchical basis response at all fine vertices for each vertex in
        rows_subset: psi[v_fine, j]. Forward pass in depth order."""
        order = np.argsort(self.depth)
        psi = np.zeros((len(self.hv), len(rows_subset)))
        colof = {r: j for j, r in enumerate(rows_subset)}
        for j, r in enumerate(rows_subset):
            psi[r, j] = 1.0
        for r in order:
            p0, p1 = self.par[r]
            if p0 < 0:
                continue
            r0, r1 = self.id2row[int(p0)], self.id2row[int(p1)]
            psi[r] += 0.5 * (psi[r0] + psi[r1])
            if r in colof:
                psi[r, colof[r]] = 1.0  # own surplus dominates; mid adds via parents
        return psi

    def ablation_table(self):
        """Value of each admitted surplus: change in held-out law-weighted
        quadratic score when the surplus is zeroed (surface-only ablation,
        no refit; metric stated in README)."""
        adm = np.where(self.admitted & (self.par[:, 0] >= 0))[0]
        if len(adm) == 0:
            return pd.DataFrame()
        psi = self.psi_matrix(list(adm))
        # response at pools: X maps fine-vertex heights -> pool values
        h2tv = {h: tv for tv, h in enumerate(self.tv2h)}
        P = np.zeros((len(self.vxy), len(adm)))
        for hr in range(len(self.hv)):
            P[h2tv[hr]] = psi[hr]
        R = self.X @ P  # pools x admitted, in height units per unit surplus
        resid = (np.log(np.clip((self.k + 0.5) / (self.n + 1), 1e-9, 1 - 1e-9)
                        / np.clip(1 - (self.k + 0.5) / (self.n + 1), 1e-9, 1 - 1e-9))
                 - self.fit_lo)
        wh = np.where(self.held & (self.pool_face >= 0),
                      1.0 / (self.v_samp + self.s2), 0.0)
        rows = []
        for j, r in enumerate(adm):
            delta = -self.surplus[r] * R[:, j] / 100.0  # logits, zeroing surplus
            dscore = 0.5 * (wh * ((resid + 0 - delta) ** 2 - resid ** 2)).sum()
            rows.append({"hier_row": int(r), "depth": int(self.depth[r]),
                         "x": float(self.hv.x[r]), "y": float(self.hv.y[r]),
                         "surplus_centi": float(self.surplus[r]),
                         "lambda_v": float(self.lam[r]),
                         "heldout_value": float(dscore)})
        return pd.DataFrame(rows).sort_values("heldout_value", ascending=False)

    # -------------------------------------------------------------- variance
    def variance_attribution(self, wala_q, wac_q):
        """Exact predictive variance of the estimate under the law, given the
        kernel: Var = sum_i k_i^2 (v_i + s2_i), decomposed by stratum."""
        kern = self.kernel(wala_q, wac_q)
        var_i = kern ** 2 * (self.v_samp + self.s2)
        sidx = np.searchsorted(STRATA_EDGES, self.n, side="right")
        tot = var_i.sum()
        return {"sd_logits": float(np.sqrt(tot)),
                "stratum_shares": [float(var_i[sidx == j].sum() / max(tot, 1e-30))
                                   for j in range(4)]}
