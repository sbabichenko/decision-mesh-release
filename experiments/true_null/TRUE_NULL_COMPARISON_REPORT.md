# True-null versus certified-null comparison

Date: 2026-07-28

## Construction

Two zero-spatial-variation benchmarks were generated on the exact Ginnie v3 coordinate and exposure design. Both use a constant logit chosen to give baseline probability 0.78.

- `clean`: binomial sampling only, no latent pool effect.
- `mix`: one latent pool effect per row from the same four size-stratified two-normal mixture law used by the deterministic histogram gate, followed by exact binomial sampling.

Ten frozen seeds were run for each benchmark through the deterministic histogram-null gate at q=0.10. No Monte Carlo was used for calibration; randomness appears only in generation of the frozen benchmark datasets.

## Results

| Benchmark | Mean admitted free coefficients | Median | Range | Runs with any admission | Mean deep faces globally | Mean fitted spatial RMS |
|---|---:|---:|---:|---:|---:|---:|
| Clean true null | 0.0 | 0 | 0--0 | 0/10 | 0.0 | 0.373 centi-logits |
| Heavy-tail true null | 2.4 | 1 | 0--11 | 6/10 | 9.4 | 8.09 centi-logits |

The clean true null shuts exactly at the root mesh in all ten runs: 81 active vertices, 128 faces, no admitted coefficients, and no depth-8 faces.

The heavy-tail true null is much harder even though the gate law matches the generator. Four of ten runs shut completely; the other six admit 1--11 free coefficients. The median fitted spatial texture remains about 7.6 centi-logits.

## Comparison with the certified v3 benchmark

Using the same histogram gate on `ginnie_certified_null_v3.csv`:

- 259 admitted free coefficients overall;
- 14 admitted coefficient centers in the certified window;
- zero admitted coefficient centers at hierarchical depth >= 8 in the window;
- 35 area-depth >= 8 faces in the window;
- prior audit: approximately 5.65 centi-logits of nonquadratic fitted texture in the window.

The true-null battery establishes that the gate can distinguish a wholly flat surface from the certified-null construction. The v3 benchmark is therefore not functioning as a generic complete-null test: most of its refinement demand comes from genuine low-order spatial variation and from representing the quadratic replacement region. Its deep-face count is a combined measure of residual false texture and geometric consequences of legitimate lower-depth fitting.

At the same time, the heavy-tail true-null battery shows that the histogram gate is not perfectly calibrated: matched heavy-tailed pool effects still create occasional false branches. Thus both conclusions hold:

1. the v3 face-count pin overstates candidate-level false discovery after B1;
2. the histogram gate still has a measurable heavy-tail false-admission problem on a genuine all-null surface.

## Recommended benchmark ledger

Use both benchmark classes:

1. **Complete null**: primary metric is admitted free coefficients and probability of any admission across seeds; secondary metrics are fitted spatial RMS and global closure footprint.
2. **Certified smooth window**: primary metrics are false admitted hierarchical surpluses and nonquadratic fitted energy after projecting out the known quadratic truth; deep-face count remains a geometry diagnostic.

The clean complete-null result is passed. The measured-heavy-tail complete-null result is not yet strong enough to claim q=0.10 frequentist FDR calibration from ten seeds, but it is far closer to the intended behavior than the raw v3 face count suggests.
