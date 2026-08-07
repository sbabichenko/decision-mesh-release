# WAC interacting-geometry EB prototype

## Model tested

The piecewise-linear mesh and quadratic child interpolation were left unchanged. For a free child vertex v, define the quadratic geometric surplus

    e_v = h_v - Q_v h.

The prototype replaced the independent EB penalty on e_v with

    r_v = e_v - rho * sum_u W_vu * sqrt(tau_v / tau_u) * e_u,

where u ranges over free vertices one depth coarser, W_vu is a triangular kernel in normalized WAC multiplied by the square root of training information, and the relation remains broad in WALA. The full sparse derivative of every r_v was included in the coupled exact solver.

The primary screened setting used a raw-WAC half-width of 0.125 percentage points and rho = 0.5.

## Seed-7 screen

On frozen topology, all tested bandwidth/strength combinations changed marginal NLL by less than 5e-5 and changed the middle-aged NLL by only a few parts in 1e-6 to 1e-5.

Allowing this parent rule during adaptive refinement was unsafe: on seed 101 it reduced the mesh from 1,042 to 206 faces and worsened overall NLL by 0.0415 and ramp NLL by 0.3455. Seed 103 also lost substantial topology.

## Matched 12-seed frozen-topology result

Compared with the quadratic-prior baseline on exactly the same fitted topologies:

- Overall marginal NLL: mean change -0.0000155; improved 7/12 seeds.
- Middle-aged NLL: mean change -0.0000026; improved 2/12 seeds.
- Ramp NLL: mean change +0.0001328; improved 4/12 seeds.
- Seasoned NLL: mean change -0.0000047; improved 10/12 seeds.

Negative change means improvement. The effect is numerically negligible and does not target the WAC ridges.

## Conclusion

This result rejects the specific rule "pool a child's quadratic surplus toward same-sign, previous-depth surpluses at nearby WAC." It does not reject interacting geometries generally. A better next formulation should introduce an explicit WAC-geometry latent effect or group mean, rather than averaging hierarchical surpluses whose signs and scales depend on refinement phase.
