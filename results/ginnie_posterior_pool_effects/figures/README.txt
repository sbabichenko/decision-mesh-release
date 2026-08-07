Corrected Ginnie Mae posterior pool-effect figures

The C++ output stores u_hat in centi-log-odds. These figures convert it to
log-odds by dividing by 100.

- Scatter: all 120,595 pools plotted; only the color scale is clipped at
  the 1st and 99th percentiles.
- Full histogram: all 120,595 posterior effects included.
- Central histogram: central 99.8% shown for readability; 121 effects
  are below the displayed range and 121 are above it.
- The top axis uses the first-order approximation 100*u percent, interpreted
  as an approximate change in termination odds for modest effects.
