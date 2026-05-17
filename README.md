# ADCS Attitude Estimation — EKF vs UKF vs Particle Filter

![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-orange?logo=mathworks)
![License](https://img.shields.io/badge/License-MIT-blue)
![Status](https://img.shields.io/badge/Status-Complete-brightgreen)
![Monte Carlo](https://img.shields.io/badge/Monte%20Carlo-50%20runs-purple)
![Filters](https://img.shields.io/badge/Filters-EKF%20%7C%20UKF%20%7C%20PF-informational)

> A rigorous, statistically validated comparison of three attitude estimation filters for a **LEO CubeSat Nadir-Pointing mission**, implemented in MATLAB with full Monte Carlo analysis.

**Author:** Hatem Adel Saber  
**Affiliation:** Space Electronics & Communications — Beni Suef University, Egypt  
**Context:** MSc Research — AI-based Attitude Determination and Control Systems (ADCS)

---

## Key Results

| Filter | Mean RMSE | Std | vs EKF | Compute |
|--------|-----------|-----|--------|---------|
| **UKF** | **0.0046°** | ±0.0004° | **−99.6%** ✅ | 6.4 s |
| EKF | 1.1233° | ±0.0056° | baseline | 0.8 s |
| PF (500p) | 12.65° | ±1.82° | −1024% ❌ | 68.4 s |

> **Statistical significance:** Wilcoxon signed-rank test, *p* = 7.56 × 10⁻¹⁰, Cohen's *d* = 281.86 (large effect) over **50 independent Monte Carlo runs**.

---

## Stress Test Results (RMSE °)

| Scenario | EKF | UKF | PF |
|----------|-----|-----|-----|
| S1 — Nominal | 1.122 | **0.004** | 11.66 |
| S2 — High IMU noise (×10) | 0.200 | **0.014** | 13.96 |
| S3 — Degraded Star Tracker (×6) | 1.107 | **0.014** | 11.02 |
| S4 — Both degraded (worst case) | 0.253 | **0.031** | 24.14 |
| S5 — Star Tracker outage | 1.591 | **1.106** | 10.24 |

> UKF outperforms EKF across **all 5 scenarios**.

---

## Mission Scenario
Mission    : LEO CubeSat — Nadir Pointing
Altitude   : 550 km circular orbit
Sim time   : 600 s  (≈ 10 min)
Time step  : 0.1 s  →  6001 steps
Init error : 12°   (realistic post-detumbling)
Seed       : rng(42) — fixed for single-run fairness

### Sensor Noise Parameters

| Sensor | Parameter | Nominal | Degraded |
|--------|-----------|---------|----------|
| IMU Gyro | σ per axis | 5×10⁻⁵ rad/s | ×10 |
| IMU Gyro | Bias drift | 1×10⁻⁶ rad/s² | ×10 |
| Star Tracker | Cross-boresight | 10 arcsec | ×6 |
| Star Tracker | Boresight/roll | 30 arcsec | ×6 |
| Star Tracker | Update rate | 1 Hz | 0.2 Hz / outage |

---

## Filter Implementations

### EKF — Extended Kalman Filter
- State: `x = [q(4×1); bias(3×1)]` — 7 states
- Linearisation via analytical Jacobians
- Joseph-form covariance update for numerical stability

### UKF — Unscented Kalman Filter
- **6-state error formulation** (avoids quaternion constraint singularity)
- Covariance `P` always 6×6 — rotation error vector + bias error
- Sigma points generated in error space, composed onto nominal quaternion
- Parameters: α = 10⁻³, β = 2, κ = 0

### PF — Particle Filter (SIR)
- 500 particles, Sequential Importance Resampling
- Systematic resampling when N_eff < N/2
- Note: performance limited by **curse of dimensionality** (7D state space)

---

## Project Structure
adcs_filter_comparison/
│
├── run_full_simulation.m              ← Entry point — runs all modules in order
│
├── shared_simulation_environment.m   ← Orbital params, true attitude, noise config
├── imu_and_star_tracker_simulation.m ← IMU gyro + Star Tracker noisy measurements
│
├── extended_kalman_filter.m          ← EKF with analytical Jacobians
├── unscented_kalman_filter.m         ← UKF with 6-state error formulation
├── particle_filter.m                 ← SIR Particle Filter (500 particles)
│
├── performance_metrics_comparison.m  ← RMSE, convergence, compute time tables
├── attitude_visualization.m          ← Euler angles, 3D trajectory, bias, 3σ envelope
├── stress_tests_and_robustness.m     ← 5 scenarios: nominal → worst case + heatmap
└── monte_carlo_statistical_analysis.m← 50-run MC, Wilcoxon test, Cohen's d, CDF plots

---

## How to Run

**Full pipeline:**
```matlab
cd adcs_filter_comparison
run_full_simulation
```

**Requirements:** MATLAB R2021a or later. No additional toolboxes required.

---

## Research Contribution

1. **Fair comparison infrastructure** — shared environment, fixed seeds, identical sensor measurements across all filters.
2. **6-state UKF formulation** — eliminates quaternion constraint singularity in 7-state covariance matrices.
3. **Monte Carlo validation** — 50 runs with Wilcoxon signed-rank testing and Cohen's *d* effect size.

---

## License

MIT License — free to use, modify, and distribute with attribution.
