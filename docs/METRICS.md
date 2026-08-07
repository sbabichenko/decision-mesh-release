# Decision Mesh: standing metric card
All values current as of the July 22 session; real data = April 2026 file.

## Predictive metrics (real data and benchmarks)
1. Held-out marginal NLL per pool (primary scoreboard).
   Fit on even half, score odd half; pool effect integrated out under a
   declared u-law. Laws: G35 (Gaussian 0.35 flat), Gstrat (Gaussian
   per-stratum total), MIX (per-stratum mixture; now standard).
   Current (MIX): default 2.978 | comp q30 2.996 | comp q10 3.034 |
   perm 3.131. Mixture improves every surface by 0.07-0.40 nats.
2. Held-out deviance per pool (secondary; mega-fragile).
   Known failure: dominated by mega pools (19 pools carry most of it),
   no credit for the pool-effect law; can invert rankings (measured:
   q10 10.8 "best" vs default 18.2 "worst" while NLL ranks reverse).
3. Mesh size (faces / vertices): parsimony indicator, reported always.

## Inference metrics (certified benchmarks only)
4. Certified-false count: faces at depth >= 8 inside the certified
   window interior (105-175 x 3.5-5.5). Current ledger:
   v3 (clean noise): default 886 | composite 2 | permutation 66
   v4 (measured physics): default 997 | composite 7 | permutation 6
   v4fs (adversarial tails): composite 2 | permutation 53 (mechanism
   open; v4fs retained as stress case, superseded as realism by the
   upcoming realistic-tails benchmark)
5. True-miss (in-window fitted-minus-truth, weighted rms, centi-logits):
   chirp: default 41.0 | composite 33.0 | permutation 30.8
   v4fs: composite 39.0 | permutation 42.8 (local excursions 59 vs 99)
6. Painting coefficient (composition-confound benchmark): share of an
   omitted-covariate field painted onto the surface: default 0.248 of
   0.35 | composite 0.183. No gate fixes this; guard covariate does.

## Parameter estimates (reported per fit, per size stratum)
7. Pool-effect law (mixture EM, exact binomial, train half), real data,
   composite q10 surface:
   n 25-128:   pi 0.280  s0^2 0.38  sL^2 2.01  total 0.84
   n 128-512:  pi 0.138  s0^2 0.28  sL^2 1.85  total 0.50
   n 512-4096: pi 0.120  s0^2 0.32  sL^2 1.42  total 0.45
   n 4096+:    pi 0.116  s0^2 0.05  sL^2 2.55  total 0.33
   Surface-dependence is a bundle diagnostic: default surface leaves
   totals 0.59/0.47/0.77/0.55; permutation surface inflates to 2.65
   (unabsorbed structure read as pool law).
8. Within-stratum residual shape: excess kurtosis 1.8 / 1.9 / 2.8 /
   20.1; P(|t|>4) = 12x / 20x / 27x / 194x Gaussian. Tails concentrate
   at mega pools.
9. Reliability lambda = sigma^2/(sigma^2+v) per stratum (drives the
   posterior kurtosis transform kurt(u_hat) = lambda^2 kurt(u)).
10. Issuer share of residual pool variance: 18.9%; multi-issuer pools
    5.4x quieter (measured, earlier session).
11. Effective replication: n_eff ceiling ~18 loans/pool (47.6M nominal
    = 4.2M effective); deep candidates run at median K_eff 2-3.

## Doctrine
12. Set-identification bounds: report default and principled fits as a
    bracket; one snapshot never completes the surface/pool split.
13. Every experiment gets a 2x2. Marginal NLL under MIX is the primary
    predictive scoreboard; certified-false counts are the primary
    inference scoreboard; deviance is a check, never a headline.
