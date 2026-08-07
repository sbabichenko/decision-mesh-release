#include "diagnostics.h"
#include "edge.h"
#include "face.h"
#include "mesh.h"
#include "observation.h"
#include "vertex.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <random>
#include <utility>
#include <vector>

namespace {
using Obs = Observation;
struct SparseSPD {
    int n = 0;
    std::vector<int> rowptr, col;
    std::vector<double> val, diag;
    // Reused CG workspace.  PM and leverage diagnostics perform many solves
    // with the same dimension; allocating four V-vectors per solve is needless.
    mutable std::vector<double> cg_r, cg_z, cg_p, cg_Ap;
    void build(int n_, const std::map<std::pair<int,int>, double>& coo) {
        n = n_;
        rowptr.assign(n + 1, 0);
        for (auto& kv : coo) rowptr[kv.first.first + 1]++;
        for (int i = 0; i < n; ++i) rowptr[i + 1] += rowptr[i];
        col.resize(coo.size()); val.resize(coo.size());
        std::vector<int> cur(rowptr.begin(), rowptr.end() - 1);
        diag.assign(n, 0.0);
        for (auto& kv : coo) {
            int r = kv.first.first, c = kv.first.second;
            int p = cur[r]++;
            col[p] = c; val[p] = kv.second;
            if (r == c) diag[r] = kv.second;
        }
    }
    void mul(const double* x, double* y) const {
        for (int r = 0; r < n; ++r) {
            double s = 0;
            for (int p = rowptr[r]; p < rowptr[r + 1]; ++p) s += val[p] * x[col[p]];
            y[r] = s;
        }
    }
    // Jacobi-preconditioned CG with caller-owned workspace.  The explicit
    // Workspace form reuses buffers across deterministic Hutchinson probes.
    int solve_ws(const double* b, double* x, std::vector<double>& vr,
                 std::vector<double>& vz, std::vector<double>& vp,
                 std::vector<double>& vAp, double tol = 1e-10,
                 int maxit = 4000, bool warm_start = false) const {
        vr.resize(n); vz.resize(n); vp.resize(n); vAp.resize(n);
        double* r = vr.data(); double* z = vz.data();
        double* p = vp.data(); double* Ap = vAp.data();
        double bn = 0;
        for (int i = 0; i < n; ++i) bn += b[i] * b[i];
        bn = std::sqrt(bn);
        if (bn == 0) { std::fill(x, x + n, 0.0); return 0; }
        if (warm_start) {
            mul(x, Ap);
            for (int i = 0; i < n; ++i) r[i] = b[i] - Ap[i];
        } else {
            std::fill(x, x + n, 0.0);
            for (int i = 0; i < n; ++i) r[i] = b[i];
        }
        double rn0 = 0;
        for (int i = 0; i < n; ++i) rn0 += r[i] * r[i];
        if (std::sqrt(rn0) < tol * bn) return 0;
        for (int i = 0; i < n; ++i) z[i] = r[i] / diag[i];
        std::copy(z, z + n, p);
        double rz = 0;
        for (int i = 0; i < n; ++i) rz += r[i] * z[i];
        for (int it = 1; it <= maxit; ++it) {
            mul(p, Ap);
            double pAp = 0;
            for (int i = 0; i < n; ++i) pAp += p[i] * Ap[i];
            double a = rz / pAp;
            double rn = 0;
            for (int i = 0; i < n; ++i) { x[i] += a * p[i]; r[i] -= a * Ap[i]; rn += r[i] * r[i]; }
            if (std::sqrt(rn) < tol * bn) return it;
            for (int i = 0; i < n; ++i) z[i] = r[i] / diag[i];
            double rz2 = 0;
            for (int i = 0; i < n; ++i) rz2 += r[i] * z[i];
            double beta = rz2 / rz; rz = rz2;
            for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
        }
        return maxit;
    }
    int solve(const double* b, double* x, double tol = 1e-10, int maxit = 4000,
              bool warm_start = false) const {
        return solve_ws(b, x, cg_r, cg_z, cg_p, cg_Ap, tol, maxit, warm_start);
    }

};

}  // namespace

// Paule-Mandel diagonal iteration: solve sum r_i^2/(1+w_i sigma) = N - tr S(sigma)
// with tr S(sigma) the smoother trace under working weights 1/(1/w+sigma) plus the
// EB surplus ridge, via Hutchinson probes on the sparse CG. Estimates the TOTAL
// extra variance sigma_a^2 + sigma_w^2 (+ zero-lag miss) with honest df accounting.
double pm_sigma_total(DecisionMesh* mesh, const std::vector<Obs>& obs,
                             const std::vector<double>& sx, const std::vector<double>& sy) {
    int N = (int)obs.size();
    // residuals and binomial info at the final surface (fixed through the iteration)
    std::vector<double> rr(N), wi(N);
    for (int i = 0; i < N; ++i) {
        double g = mesh->predict_index(i);
        double ph = std::clamp(1 / (1 + std::exp(-(obs[i].offset + g / 100.0))), 1e-9, 1 - 1e-9);
        wi[i] = obs[i].exposure * ph * (1 - ph);
        rr[i] = (obs[i].events - obs[i].exposure * ph) / std::sqrt(wi[i]);
    }
    // design triplets on the final mesh
    std::vector<Vertex*> verts;
    for (Vertex* v : mesh->vertices) if (v->active) verts.push_back(v);
    std::sort(verts.begin(), verts.end(), [](Vertex* a, Vertex* b) {
        return a->x != b->x ? a->x < b->x : a->y < b->y; });
    int V = (int)verts.size();
    if (V == 0) return 0.0;
    std::map<Vertex*, int> vid;
    for (int i = 0; i < V; ++i) vid[verts[i]] = i;
    std::vector<std::array<int, 3>> tv(N, {{-1, -1, -1}});
    std::vector<std::array<float, 3>> tw(N);
    for (Face* f : mesh->active_faces) {
        if (f->coords_indices.empty()) continue;
        Vertex* v0 = f->vertices[0]; Vertex* v1 = f->vertices[1]; Vertex* v2 = f->vertices[2];
        auto i0 = vid.find(v0), i1 = vid.find(v1), i2 = vid.find(v2);
        if (i0 == vid.end() || i1 == vid.end() || i2 == vid.end()) continue;
        double ax = v1->x - v0->x, ay = v1->y - v0->y;
        double bx = v2->x - v0->x, by = v2->y - v0->y;
        double d00 = ax * ax + ay * ay, d11 = bx * bx + by * by;
        if (d00 < 1e-14 || d11 < 1e-14) continue;
        for (int i : f->coords_indices) {
            double px = sx[i] - v0->x, py = sy[i] - v0->y;
            double u = (px * ax + py * ay) / d00, w = (px * bx + py * by) / d11;
            tv[i] = {i0->second, i1->second, i2->second};
            tw[i] = {static_cast<float>(1.0 - u - w),
                     static_cast<float>(u), static_cast<float>(w)};
        }
    }
    // ridge stencils (fixed across sigma)
    std::vector<std::array<int, 3>> rid; std::vector<double> rlam;
    for (Vertex* v : verts) {
        if (!v->parent_edge) continue;
        auto it = mesh->tau_sq.find(v->depth);
        if (it == mesh->tau_sq.end() || it->second <= 0) continue;
        auto a0 = vid.find(v->parent_edge->vertex0), a1 = vid.find(v->parent_edge->vertex1);
        if (a0 == vid.end() || a1 == vid.end()) continue;
        rid.push_back({vid[v], a0->second, a1->second});
        rlam.push_back(1.0 / it->second);
    }
    // Build the sparse matrix pattern once.  PM changes only the observation
    // weights as sigma moves, so rebuilding a std::map and CSR structure at every
    // root evaluation is pure overhead.
    std::map<std::pair<int,int>, double> pattern;
    for (int i = 0; i < N; ++i) if (tv[i][0] >= 0)
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
            pattern[{tv[i][a], tv[i][b]}] = 0.0;
    for (const auto& r : rid)
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
            pattern[{r[a], r[b]}] = 0.0;
    for (int i = 0; i < V; ++i) pattern[{i, i}] = 0.0;

    std::map<std::pair<int,int>, int> pos;
    int kk = 0;
    for (const auto& kv : pattern) pos[kv.first] = kk++;
    SparseSPD A;
    A.build(V, pattern);
    std::vector<std::array<int, 9>> obs_pos(N);
    for (int i = 0; i < N; ++i) {
        obs_pos[i].fill(-1);
        if (tv[i][0] < 0) continue;
        int q = 0;
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
            obs_pos[i][q++] = pos[{tv[i][a], tv[i][b]}];
    }
    std::vector<std::array<int, 9>> ridge_pos(rid.size());
    for (size_t k = 0; k < rid.size(); ++k) {
        int q = 0;
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
            ridge_pos[k][q++] = pos[{rid[k][a], rid[k][b]}];
    }
    std::vector<int> diag_pos(V);
    for (int i = 0; i < V; ++i) diag_pos[i] = pos[{i, i}];
    constexpr int PM_PROBES = 6;
    std::vector<double> wt(N);
    // Precompute the exact deterministic probe stream once.  Six int8 signs
    // per observation cost little and remove PRNG work from every PM evaluation.
    std::vector<signed char> probe_sign((size_t)PM_PROBES * N);
    { std::mt19937 prng(999); for (size_t j = 0; j < probe_sign.size(); ++j)
          probe_sign[j] = (prng() & 1) ? 1 : -1; }
    // One warm-start vector per deterministic Hutchinson probe.  The PM
    // bisection changes sigma smoothly, so these are very close to the next
    // solutions and substantially reduce CG iterations.
    std::vector<double> probe_x((size_t)PM_PROBES * V, 0.0);
    std::vector<unsigned char> probe_seen(PM_PROBES, 0);
    static const double co[3] = {1.0, -0.5, -0.5};

    auto trace_S = [&](double sig, int nprobe) {
        std::fill(A.val.begin(), A.val.end(), 0.0);
        for (int i = 0; i < N; ++i) {
            if (tv[i][0] < 0) { wt[i] = 0; continue; }
            wt[i] = 1.0 / (1.0 / wi[i] + sig) / 1e4;   // match ridge units
            int q = 0;
            for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
                A.val[obs_pos[i][q++]] += wt[i] * tw[i][a] * tw[i][b];
        }
        for (size_t k = 0; k < rid.size(); ++k) {
            int q = 0;
            for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
                A.val[ridge_pos[k][q++]] += rlam[k] * co[a] * co[b];
        }
        for (int i = 0; i < V; ++i) A.val[diag_pos[i]] += 1e-8;
        for (int i = 0; i < V; ++i) A.diag[i] = A.val[diag_pos[i]];

        double tr = 0;
        // Deliberately serial: fixed probe order and accumulation order make
        // the diagnostic bitwise reproducible on a fixed toolchain/CPU.
        for (int p = 0; p < nprobe; ++p) {
            std::vector<double> rhs(V, 0.0), wr, wz, wp, wAp;
            const signed char* ev = probe_sign.data() + (size_t)p * N;
            for (int i = 0; i < N; ++i) {
                if (tv[i][0] < 0) continue;
                double we = wt[i] * (double)ev[i];
                for (int a = 0; a < 3; ++a) rhs[tv[i][a]] += we * tw[i][a];
            }
            double* xp = probe_x.data() + (size_t)p * V;
            A.solve_ws(rhs.data(), xp, wr, wz, wp, wAp, 1e-8, 2000,
                       probe_seen[p] != 0);
            probe_seen[p] = 1;
            double acc = 0;
            for (int i = 0; i < N; ++i) {
                if (tv[i][0] < 0) continue;
                double Se = 0;
                for (int a = 0; a < 3; ++a) Se += tw[i][a] * xp[tv[i][a]];
                acc += (double)ev[i] * Se;
            }
            tr += acc;
        }
        return tr / nprobe;
    };
    auto srr_at = [&](double sig) {
        double srr = 0;
        for (int i = 0; i < N; ++i) srr += rr[i] * rr[i] / (1.0 + wi[i] * sig);
        return srr;
    };
    // Given a fixed trace, solve the scalar PM equation exactly (up to a cheap
    // bisection tolerance). This avoids recomputing six Hutchinson solves at
    // every root midpoint. The smoother trace changes slowly with sigma, so a
    // few deterministic fixed-point updates suffice.
    auto solve_given_trace = [&](double tr) {
        const double target = (double)N - tr;
        if (srr_at(0.0) <= target) return 0.0;
        double lo = 0.0, hi = 1.0;
        while (srr_at(hi) > target && hi < 8.0) hi *= 2.0;
        for (int it = 0; it < 36; ++it) {
            double mid = 0.5 * (lo + hi);
            if (srr_at(mid) > target) lo = mid; else hi = mid;
        }
        return 0.5 * (lo + hi);
    };

    // Start at the binomial-scale smoother trace, then update trace and sigma.
    // Five evaluations cap PM at 30 sparse solves while preserving the scalar root accuracy.
    double tr = trace_S(0.0, PM_PROBES);
    double sol = solve_given_trace(tr);
    int neval = 1;
    for (; neval < 5; ++neval) {
        double tr_new = trace_S(sol, PM_PROBES);
        double sol_new = solve_given_trace(tr_new);
        tr = tr_new;
        if (neval >= 3 && std::abs(sol_new - sol) <= 1e-8 * std::max(1.0, sol)) {
            sol = sol_new;
            ++neval;
            break;
        }
        sol = sol_new;
    }
    std::fprintf(stderr, "[trS] %.2f at sigma %.4f (%d trace evals)\n", tr, sol, neval);
    return sol;
}

// Analytic leverage-corrected sigma_u^2 on the final mesh.
// Returns corrected estimate; also fills D_sig/D_naive ratio and C_noise/T.
double leverage_correction(DecisionMesh* mesh, const std::vector<Obs>& obs,
                                  const std::vector<double>& sx, const std::vector<double>& sy,
                                  int npool, double naive_sig, double& dsig_ratio,
                                  double& cnoise_ratio, int leverage_probes) {
    // index active vertices
    std::vector<Vertex*> verts;
    for (Vertex* v : mesh->vertices) if (v->active) verts.push_back(v);
    std::sort(verts.begin(), verts.end(), [](Vertex* a, Vertex* b) {
        return a->x != b->x ? a->x < b->x : a->y < b->y; });
    int V = (int)verts.size();
    if (V == 0) { dsig_ratio = 1; cnoise_ratio = 0; return naive_sig; }
    std::map<Vertex*, int> vid;
    for (int i = 0; i < V; ++i) vid[verts[i]] = i;

    // design triplets from active faces: obs index -> (3 vids, 3 barycentric weights)
    int N = (int)obs.size();
    // per-obs binomial information at the FINAL surface (matches naive stage units)
    std::vector<float> infoN(N);
    for (int i = 0; i < N; ++i) {
        double g = mesh->predict_index(i);
        double ph = 1 / (1 + std::exp(-(obs[i].offset + g / 100.0)));
        infoN[i] = obs[i].exposure * ph * (1 - ph);
    }
    std::vector<std::array<int, 3>> tv(N, {{-1, -1, -1}});
    std::vector<std::array<float, 3>> tw(N);
    for (Face* f : mesh->active_faces) {
        if (f->coords_indices.empty()) continue;
        Vertex* v0 = f->vertices[0]; Vertex* v1 = f->vertices[1]; Vertex* v2 = f->vertices[2];
        auto i0 = vid.find(v0), i1 = vid.find(v1), i2 = vid.find(v2);
        if (i0 == vid.end() || i1 == vid.end() || i2 == vid.end()) continue;
        double ax = v1->x - v0->x, ay = v1->y - v0->y;
        double bx = v2->x - v0->x, by = v2->y - v0->y;
        double d00 = ax * ax + ay * ay, d11 = bx * bx + by * by;
        if (d00 < 1e-14 || d11 < 1e-14) continue;
        for (int i : f->coords_indices) {
            double px = sx[i] - v0->x, py = sy[i] - v0->y;
            double u = (px * ax + py * ay) / d00;
            double w = (px * bx + py * by) / d11;
            tv[i] = {i0->second, i1->second, i2->second};
            tw[i] = {static_cast<float>(1.0 - u - w),
                     static_cast<float>(u), static_cast<float>(w)};
        }
    }

    // A = Phi' W~ Phi assembled sparse (COO -> CSR); Bm no longer materialized
    // (noise channel verified negligible and gated off).
    const std::vector<double>& wfit = mesh->weights;
    std::map<std::pair<int,int>, double> Acoo;
    for (int i = 0; i < N; ++i) {
        if (tv[i][0] < 0) continue;
        double wf = wfit[i];
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
            Acoo[{tv[i][a], tv[i][b]}] += wf * tw[i][a] * tw[i][b];
    }
    int ridge_n = 0; double tau_mean = 0;
    for (Vertex* v : verts) {
        if (!v->parent_edge) continue;
        auto it = mesh->tau_sq.find(v->depth);
        if (it == mesh->tau_sq.end() || it->second <= 0) continue;
        double lam = 1.0 / it->second;
        Vertex* p0 = v->parent_edge->vertex0; Vertex* p1 = v->parent_edge->vertex1;
        auto a0 = vid.find(p0), a1 = vid.find(p1);
        if (a0 == vid.end() || a1 == vid.end()) continue;
        ridge_n++; tau_mean += it->second;
        int id[3] = {vid[v], a0->second, a1->second};
        double co[3] = {1.0, -0.5, -0.5};
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
            Acoo[{id[a], id[b]}] += lam * co[a] * co[b];
    }
    for (int i = 0; i < V; ++i) Acoo[{i, i}] += 1e-8;
    std::fprintf(stderr, "[corr diag] V=%d ridge_verts=%d mean_tau_sq=%.4f nnz=%zu\n",
                 V, ridge_n, ridge_n ? tau_mean / ridge_n : 0.0, Acoo.size());
    SparseSPD A;
    A.build(V, Acoo);


    // Pool loading vectors are not needed by the diagonal-leverage approximation.
    // The previous dense npool x V allocation could dominate memory (hundreds of MB).

    // Single factorization already done. Compute A^{-1} once by solving against
    // the V basis columns (O(V^3), one time), then every per-point self-leverage
    // S_ii = W~_i * phi_i' Ainv phi_i is an O(1) 3x3 form. This replaces the
    // O(npool * V^2) per-pool solves that dominated runtime.
    // Self-leverage diag(S), S = Phi A^{-1} Phi' W~. Forming A^{-1} fully is O(V^3)
    // and became a third of runtime. Instead estimate diag(S) by Hutchinson:
    //   diag(S)_v = E[ z_v (S z)_v ],  z Rademacher.
    // Each probe costs one solve (O(V^2) via the existing factorization) plus two
    // sparse Phi products (O(N)). A handful of probes suffices because per-vertex
    // leverage is smooth and we only need pool-mean accuracy.
    // S z at vertex space: g = A^{-1} (Phi' W~ Phi z_field)? Careful: S acts on the
    // data-space; we need diag in DATA space (per observation i). So probe in data
    // space: for Rademacher e (length N), (S e)_i = phi_i' A^{-1} Phi' W~ e.
    //   rhs_v = sum_j w~_j phi_j(v) e_j   (sparse accumulate, O(N))
    //   x = A^{-1} rhs                     (one solve)
    //   (S e)_i = sum_a tw[i][a] x[tv_i_a]
    // diag(S)_i = E[ e_i (S e)_i ].
    const int NPROBE = std::clamp(leverage_probes, 1, 64);
    std::vector<float> Sii_acc(N, 0.0f);
    std::mt19937 prng(12345);
    std::vector<double> rhs(V), xsol(V);
    for (int pr = 0; pr < NPROBE; ++pr) {
        std::fill(rhs.begin(), rhs.end(), 0.0);
        // Generate probe signs on demand; no N-length probe vector is retained.
        // Save the PRNG state so the same signs can be replayed after the solve.
        const std::mt19937 probe_start = prng;
        for (int i = 0; i < N; ++i) {
            const double ei = (prng() & 1) ? 1.0 : -1.0;
            if (tv[i][0] < 0) continue;
            double we = wfit[i] * ei;
            for (int a = 0; a < 3; ++a) rhs[tv[i][a]] += we * tw[i][a];
        }
        A.solve(rhs.data(), xsol.data());
        prng = probe_start;
        for (int i = 0; i < N; ++i) {
            const double ei = (prng() & 1) ? 1.0 : -1.0;
            if (tv[i][0] < 0) continue;
            double Se = 0;
            for (int a = 0; a < 3; ++a) Se += (double)tw[i][a] * xsol[tv[i][a]];
            Sii_acc[i] += (float)(ei * Se);
        }
    }
    std::vector<double> aP(npool, 0.0), bP(npool, 0.0);
    double hsum = 0; int hcnt2 = 0;
    for (int i = 0; i < N; ++i) {
        if (tv[i][0] < 0) continue;
        int P = obs[i].pool_id;
        double Sii = Sii_acc[i] / NPROBE;
        double s = std::sqrt(infoN[i]), one_h = 1.0 - Sii;
        aP[P] += s * one_h; bP[P] += infoN[i] * one_h * one_h;
        hsum += Sii; hcnt2++;
    }
    std::fprintf(stderr, "[corr diag] mean S_ii = %.4f (n=%d, %d probes)\n",
                 hsum/std::max(hcnt2,1), hcnt2, NPROBE);
    double D_sig = 0, D_naive = 0;
    {
        std::vector<double> sA(npool, 0.0), sB(npool, 0.0);
        for (int i = 0; i < N; ++i) {
            if (tv[i][0] < 0) continue;
            double s = std::sqrt(infoN[i]);
            sA[obs[i].pool_id] += s; sB[obs[i].pool_id] += infoN[i];
        }
        for (int P = 0; P < npool; ++P) {
            D_sig += aP[P] * aP[P] - bP[P];
            D_naive += sA[P] * sA[P] - sB[P];
        }
    }

    double C_noise = 0.0;  // noise channel negligible in DL regime (verified)


    // T (naive pairwise numerator) on final residuals
    std::vector<double> Sr(npool, 0), Sr2(npool, 0);
    for (int i = 0; i < N; ++i) {
        double g = mesh->predict_index(i);
        double ph = 1 / (1 + std::exp(-(obs[i].offset + g / 100.0)));
        double inf = obs[i].exposure * ph * (1 - ph);
        double r = (obs[i].events - obs[i].exposure * ph) / std::sqrt(inf);
        Sr[obs[i].pool_id] += r; Sr2[obs[i].pool_id] += r * r;
    }
    double T = 0;
    for (int P = 0; P < npool; ++P) T += Sr[P] * Sr[P] - Sr2[P];

    dsig_ratio = D_naive > 0 ? D_sig / D_naive : 1.0;
    cnoise_ratio = T != 0 ? C_noise / T : 0.0;
    return (T - C_noise) / std::max(D_sig, 1e-9);
}
