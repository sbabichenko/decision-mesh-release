# A pooled-permutation admission gate for Decision Mesh

Working memo — independent evaluation, August 2026.
All numbers below are from simulation on the real `data/ginnie_design.csv` exposure
geometry, using the release's own measured pool-effect law (`core/auditing/law_shared.json`).

---

## 1. The problem

The shipped candidate gate prices its null variance from a *model*:

```
sd² = (1/X_v + σ²_gate · X⁽²⁾_v / X_v² + s_p/4) · φ̂
```

Every ingredient — `σ²_gate = 0.35`, the per-stratum variance table, the kurtosis
table — is estimated from the same data the gate is judging. The mega stratum of
the fitted law has kurtosis 21.

Measured null size against a flat truth (nominal 0.05):

| scenario | shipped gate |
|---|---|
| clean binomial, p = 0.78 | 0.000 |
| clean binomial, p = 0.01 | 0.019 |
| **Gaussian pool effect, p = 0.01** | **0.157** |
| heavy-tail law, p = 0.01 | 0.142 |
| issuer-clustered | 0.188 |

Two things matter here. The gate is simultaneously *over*-conservative when there is
no pool heterogeneity and 3× *anti*-conservative when there is. And it fails at 0.157
under a purely **Gaussian** pool law — so this is a wrong-variance problem, not a
heavy-tails problem, and no amount of kurtosis pricing addresses it.

End-to-end confirmation on the release binary: a matched heavy-tail complete null
(true surface exactly constant) admits ~450 false vertices per run at p = 0.01,
10/10 seeds, with 0.70 log-odds of manufactured surface.

### Diagnosis

Two channels, one root cause.

1. **K_eff ≈ 1.** Admitted vertices have median effective pool count 1.3 — the
   statistic is a single draw from the pool-effect law. No moment correction can
   help; there is no averaging for a CLT to act on. This is why the one-law gate
   (which *does* price kurtosis) performed no better than flat composite.
2. **The empirical null is fitted on the wrong scale.** Central matching sees a
   score family whose centre is a mass of near-zero-information candidates, fits
   σ̂₀ ≈ 1, and prices a tail whose true sd is 1.6–7.

Root cause: the null is modelled rather than measured from the candidate's own pools.

---

## 2. The approach

Three components. The third is what makes it work.

### 2.1 Permutation statistic

Under H₀ (no spatial signal in the star) the pairs (residual, variance) are
exchangeable across geometric slots. With `a` the (geometric) hat-column coefficients
and `gᵢ = wᵢeᵢ` the pool score contributions, the statistic is

```
N = Σᵢ aᵢ g_π(i)
```

with `g` a fixed multiset assigned at random to slots. No pool law is required.

The self-normalising property is the point: for a star dominated by one pool the
studentised version tends to ±1 regardless of the excursion size, so a single-pool
star can never be significant. The minimum-effective-pools rule becomes a
consequence rather than a hand-set floor.

### 2.2 Exact moments, not Monte Carlo

The randomisation moments are closed form (Pitman). With `A_r`, `B_r` the centred
power sums of `a` and `g`:

```
E[N]  = m ā ḡ
Var[N] = A₂B₂ / (m−1)
μ₃[N]  = A₃B₃ / ((m−1)(m−2))
```

For tail probabilities, coarsening `a` to two levels reduces the permutation law to a
fixed-size **subset sum**, which admits a double saddlepoint. Validated against
complete enumeration of C(m,r) assignments: median ratio 1.03 at exact p ≤ 0.10;
median relative error 1.6% overall.

Saddlepoint beats Cornish–Fisher decisively in the tail, and the gap widens where it
matters (ground truth = enumeration):

| exact tail | Gaussian | Cornish–Fisher | saddlepoint |
|---|---|---|---|
| 0.050 | 0.836× | 0.962× | **0.973×** |
| 0.010 | 1.329× | 1.539× | **0.998×** |
| 0.002 | 2.022× | 2.659× | **0.978×** |

### 2.3 Cross-star pooling — the essential step

A within-star permutation cannot produce a p-value below 1/C(m,r). Median across real
stars is 2.0e-3, and **35% of stars cannot reach p = 0.01 at all**. Ignoring that
boundary makes the continuous saddlepoint 17–34× anti-conservative below 1e-3.

The fix is not to floor but to pool. Standardise each star by its **own** exact
moments (Cornish–Fisher, using E, V, μ₃ above), then pool the *residual shape* across
the family. That shape is near-homogeneous:

| stratum | P(\|z\|>2) | P(\|z\|>3) |
|---|---|---|
| K_eff < 2 | 0.0281 | 0.0008 |
| K_eff 2–6 | 0.0311 | 0.0012 |
| K_eff ≥ 6 | 0.0325 | 0.0006 |

Build the null reference from permutation replicates across all stars: M stars × B
draws, resolution 1/(M·B), banded by star size m.

---

## 3. Results

### Null size, all deterministic rules (nominal 0.05 / 0.01)

| scenario | shipped z | perm saddlepoint | donor ×4 |
|---|---|---|---|
| clean binomial p = 0.78 | 0.000 / 0.000 | 0.049 / 0.020 | 0.043 / 0.007 |
| clean binomial p = 0.01 | 0.019 / 0.005 | 0.064 / 0.022 | 0.045 / 0.007 |
| Gaussian pool p = 0.01 | 0.157 / 0.083 | 0.045 / 0.017 | 0.035 / 0.008 |
| heavy-tail law p = 0.78 | 0.095 / 0.053 | 0.069 / 0.026 | 0.044 / 0.007 |
| heavy-tail law p = 0.01 | 0.142 / 0.092 | 0.055 / 0.020 | 0.066 / 0.033 |
| heavy, mega stars | 0.082 / 0.062 | 0.057 / 0.023 | 0.068 / 0.013 |
| issuer-clustered | 0.188 / 0.126 | 0.045 / 0.021 | 0.054 / 0.012 |

The permutation rules hold 0.035–0.069 across every scenario while the shipped gate
ranges 0.000–0.188. Notably the issuer-clustered scenario (an exchangeability
violation) does not break the permutation gate.

### Deep-tail calibration with the pooled null

| nominal | pooled null | theoretical N(0,1) |
|---|---|---|
| 0.05 | **0.0500** | 0.0368 |
| 0.01 | 0.0107 | 0.0051 |
| 0.001 | 0.00129 | 0.00071 |
| 0.0001 | 0.00007 | 0.00007 |

Smallest achievable p: **1.6e-5** vs the per-star floor of 2.0e-3 — a 125×
improvement, extensible by increasing B.

Contamination check (pooled null built from a family that is 15% signal): calibration
on the true-null members holds at 0.052 / 0.0090 / 0.00134.

### Multiplicity and the q dial

| rule | q=0.05 | 0.10 | 0.20 | 0.30 | 0.40 | mean \|FDP−q\| |
|---|---|---|---|---|---|---|
| BH | 0.13/0.01 | 0.12/0.02 | 0.19/0.06 | 0.25/0.10 | 0.34/0.15 | 0.044 |
| Storey-BH | 0.13/0.01 | 0.11/0.02 | 0.20/0.07 | 0.26/0.10 | 0.36/0.17 | 0.035 |
| **lfdr, theoretical null** | 0.00/0.01 | 0.10/0.03 | 0.14/0.10 | 0.28/0.17 | 0.39/0.25 | **0.028** |

(cells are FDP / power)

BH becomes *viable* only after cross-star pooling — before it, BH thresholds q/M lie
below the minimum attainable p for 75–80% of candidates. It remains less powerful
than lfdr. Using the **theoretical** null rather than an empirical one roughly halves
the calibration error, because permutation p-values are uniform by construction.

### Power, size-adjusted, by replication

| band | power at δ=2 |
|---|---|
| K_eff < 2 | 0.055 |
| K_eff ≥ 2 | **0.433** |

Power is concentrated almost entirely in replicated stars — the under-replicated band
contributes risk without evidence.

### Cost per candidate

| | µs |
|---|---|
| shipped composite z | 6.1 |
| CR sandwich t | 14.3 |
| perm saddlepoint (scipy) | 897 |
| perm saddlepoint (hand Newton) | 529 |
| donor saddlepoint (1-D) | 728 |

The hand-rolled Newton version agrees with scipy to 1e-14. In C++ with warm starts
this should land within ~10× of the shipped statistic. Exact moments alone are O(m)
and effectively free, so they can serve as a pre-screen.

---

## 4. Negative results worth recording

- **Pooling the variance fails.** A limma-style moderated statistic (shrink each
  star's sandwich variance toward a pooled prior) gives size 0.195 at ν₀=5 and 0.273
  at ν₀=20 — worse than the shipped gate. A star's variance is determined largely by
  whether a lump landed in it, so shrinking toward the median under-states variance
  for exactly the stars that generate false admissions.
- **Donor bank size is a validity knob, not a power knob.** Augmenting the reference
  with matched donors from elsewhere doubles power at ×4, but degrades sharply beyond:
  size at nominal 0.001 is 0.000 (×4), 0.001 (×25), **0.078** (×100). Past ×4 the
  reference becomes the global law and reinstates the estimated-law failure.
- **Monotone (isotonic) lfdr** does not improve q-calibration; it costs power at small q.
- **A K_eff floor helps but saturates.** Raising `DMESH_MIN_CURRENT_KEFF` from 0.5 to
  32 cuts null admissions ~90× (445 → 5) but never reaches zero. Separately, the guard
  is checked at candidacy only: vertices admitted at a floor of 8 end up with median
  final K_eff ≈ 2, because descendants subdivide the star afterwards.

---

## 5. Proposed rule

1. **Screen** on exact permutation moments and the CR sandwich (O(m), ~15 µs).
2. **Test** with the two-level permutation saddlepoint (~530 µs, deterministic).
3. **Calibrate** by lookup against a cross-star pooled null built from permutation
   replicates, banded by star size m.
4. **Admit** by lfdr prefix against the **theoretical** N(0,1) null; run central
   matching only as an audit that should return σ̂₀ ≈ 1.

Independent of the above, two changes are free and worth making immediately:
raise `DMESH_MIN_CURRENT_KEFF` from 0.5 to ~4 (halves null admissions for +0.003 NLL),
and re-check K_eff at commit rather than only at candidacy.

---

## 6. Limitations

- **Fixed topology, single round.** Everything here is a per-candidate null. The
  sequential channel — descendants of a false admission getting repeated shots at the
  same realised anomaly — is untouched and still needs branchwise error spending.
- **Emulated geometry.** Stars are grid cells with a 1-D hat column, not NVB
  hierarchical columns with EB shrinkage inside the loop. The real columns are
  continuous barycentric values; the two-level coarsening is a design choice whose
  power cost was measured but not optimised.
- **Power ceiling.** At controlled FDP, power tops out around 0.25–0.43. That is the
  single-cross-section identification bound, not a defect of the machinery. A panel
  would change it.
- **Band resolution.** The deepest reachable p is 1/(M·B) — a compute knob. For a
  24,404-candidate family, B = 100 gives ~4e-7.
- **Stratification is coarse.** Bands by m are needed (small stars run conservative,
  m 10–20 runs mildly hot); finer bands or a smooth kernel over m would help.

## 7. Corrections made during this work

Two bugs were found and fixed in the evaluation code itself; both had produced
plausible-looking but invalid numbers:

- The two-sided p-value was assembled as `P(N≥t) + P(N≥2E−t)` instead of
  `2·min(P(N≥t), 1−P(N≥t))`. The broken form returns p ≈ 1 for 88% of calls and
  silently falls back to a Gaussian for the rest; the two errors cancelled to
  something that resembled nominal calibration.
- An lfdr comparison was run on unsigned z, which breaks the two-groups model.

Both were caught by validating against exact enumeration rather than against another
approximation. Any further work on this should keep enumeration-based unit tests as
the ground truth — for m ≤ 20 they are milliseconds, and they are the only check that
does not share assumptions with the thing being tested.