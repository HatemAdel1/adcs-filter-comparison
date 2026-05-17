# ADCS Filter Comparison Project
## EKF vs UKF vs Particle Filter — Attitude Estimation

**Author:** Hatem Adel Saber  
**MSc:** Space Electronics & Communications — Beni Suef University  
**Topic:** AI-based ADCS using Star Tracker & IMU Sensor Fusion

---

## Project Structure

```
adcs_filter_comparison/
├── RUN_ALL.m               ← Start here — runs everything in order
├── M0_shared_environment.m ← Orbital parameters, true attitude, noise config
├── M1_sensor_simulation.m  ← IMU & Star Tracker noisy measurements
├── M2a_EKF.m               ← Extended Kalman Filter (with Jacobians)
├── M2b_UKF.m               ← Unscented Kalman Filter (sigma points)
├── M2c_PF.m                ← Particle Filter (SIR, 500 particles)
├── M3_comparison.m         ← RMSE, convergence, compute time comparison
├── M4_visualization.m      ← Euler angles, 3D trajectory, bias plots
└── M5_stress_tests.m       ← 5 scenarios: nominal → worst case
```

---

## Scenario

| Parameter       | Value                         |
|----------------|-------------------------------|
| Mission        | LEO CubeSat — Nadir Pointing  |
| Altitude       | 550 km circular orbit         |
| Sim time       | 600 s (10 min)                |
| Time step      | 0.1 s                         |
| Initial error  | 12° (post-detumbling)         |
| Random seed    | 42 (fixed — fair comparison)  |

---

## Noise Parameters

### IMU Gyroscope (Nominal — S1)
- σ_gyro = 5×10⁻⁵ rad/s per axis
- Bias drift = 1×10⁻⁶ rad/s²

### Star Tracker (Nominal — S1)
- Cross-boresight: 10 arcsec (48 µrad)
- Boresight/roll: 30 arcsec (145 µrad)
- Update rate: 1 Hz

---

## Stress Test Scenarios

| Scenario | IMU noise | Star Tracker | Purpose             |
|----------|-----------|--------------|---------------------|
| S1       | Nominal   | Nominal      | Baseline accuracy   |
| S2       | ×10       | Nominal      | Gyro degradation    |
| S3       | Nominal   | ×6           | Star tracker loss   |
| S4       | ×10       | ×6           | Worst case          |
| S5       | Medium    | 10s outage/min | Filter robustness |

---

## How to Run

1. Open MATLAB
2. `cd` into the project folder
3. Run: `RUN_ALL`

Or run modules individually in order: M0 → M1 → M2a → M2b → M2c → M3 → M4 → M5

---

## State Vector

```
x = [q1; q2; q3; q4; bx; by; bz]   (7×1)
    quaternion (vector+scalar) + gyro bias
```

## Key Outputs

- `EKF.rmse_deg` — EKF steady-state RMSE [deg]
- `UKF.rmse_deg` — UKF steady-state RMSE [deg]  
- `PF.rmse_deg`  — PF steady-state RMSE [deg]
- `STRESS.rmse_table` — 5×3 matrix (scenarios × filters)
