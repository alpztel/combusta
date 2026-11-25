# Benchmarks (local workspace)

This directory is a local workspace for open reference data used to compare simulations against experiments, for example the spray penetration measurements of the Sandia [Engine Combustion Network](https://ecn.sandia.gov/) (ECN), such as Spray A (n-dodecane, 0.090 mm nozzle) and Spray D (0.186 mm nozzle).

The repository ships **no** reference data. ECN files are distributed by Sandia under their own terms; download them from the ECN website and keep them here locally. Data files in this directory (`*.csv`, `*.txt`, `*.dat`) are ignored by Git so that they are not committed by accident.

## Comparing a CFD spray with ECN data

```python
from combusta import compute_spray_validation_metrics, load_ecn_spray_data, parse_spray_csv

ecn = load_ecn_spray_data("benchmarks/ECN_SprayA_900K_22.8kgm3.txt", condition_label="ECN Spray A")
cfd = parse_spray_csv("path/to/forte_spray.csv")

metrics = compute_spray_validation_metrics(
    cfd.vapor_penetration_m,
    ecn.vapor_penetration_m,
    cfd_time_s=cfd.time_s,
    ref_time_s=ecn.time_s,
)
print(f"RMSE = {metrics.rmse_m * 1e3:.2f} mm, MAPE = {metrics.mape_pct:.2f} %, "
      f"L2 = {metrics.l2_relative_error:.4f}")
```

The file name above is an example; use the files you downloaded. `load_ecn_spray_data` accepts comma- or whitespace-separated tables with `#` or `%` comment headers and converts time and length to SI units. When the ECN file contains standard deviations, they are available as `liquid_std_m` and `vapor_std_m` and can enter an ASME V&V 20 comparison (`combusta.validate`) as the experimental uncertainty.
