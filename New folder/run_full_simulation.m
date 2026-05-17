% =========================================================================
% run_full_simulation.m
% ADCS Filter Comparison Project — Master Script
%
% Runs all modules in sequence:
%   shared_simulation_environment
%   → imu_and_star_tracker_simulation
%   → extended_kalman_filter
%   → unscented_kalman_filter
%   → particle_filter
%   → performance_metrics_comparison
%   → attitude_visualization
%   → stress_tests_and_robustness
%
% Author  : Hatem Adel Saber
% MSc     : Space Electronics & Communications — Beni Suef University
% =========================================================================

clear; clc; close all;

fprintf('============================================================\n');
fprintf('  ADCS FILTER COMPARISON PROJECT\n');
fprintf('  EKF vs UKF vs Particle Filter — LEO CubeSat\n');
fprintf('  Author: Hatem Adel Saber\n');
fprintf('============================================================\n\n');

%% MODULE 0 — Shared environment
fprintf('>>> Running Module 0: Shared Environment\n');
run('shared_simulation_environment.m');
fprintf('>>> Shared environment complete.\n\n');

%% SENSOR SIMULATION
fprintf('>>> Running: IMU & Star Tracker Simulation\n');
run('imu_and_star_tracker_simulation.m');
fprintf('>>> Sensor simulation complete.\n\n');

%% EXTENDED KALMAN FILTER
fprintf('>>> Running: Extended Kalman Filter\n');
run('extended_kalman_filter.m');
fprintf('>>> EKF complete.\n\n');

%% UNSCENTED KALMAN FILTER
fprintf('>>> Running: Unscented Kalman Filter\n');
run('unscented_kalman_filter.m');
fprintf('>>> UKF complete.\n\n');

%% PARTICLE FILTER
fprintf('>>> Running: Particle Filter\n');
run('particle_filter.m');
fprintf('>>> Particle Filter complete.\n\n');

%% PERFORMANCE COMPARISON
fprintf('>>> Running: Performance Metrics & Comparison\n');
run('performance_metrics_comparison.m');
fprintf('>>> Comparison complete.\n\n');

%% VISUALIZATION
fprintf('>>> Running: Attitude Visualization\n');
run('attitude_visualization.m');
fprintf('>>> Visualization complete.\n\n');

%% STRESS TESTS
fprintf('>>> Running: Stress Tests & Robustness Analysis\n');
run('stress_tests_and_robustness.m');
fprintf('>>> Stress tests complete.\n\n');

%% FINAL SUMMARY
fprintf('============================================================\n');
fprintf('  ALL MODULES COMPLETE\n');
fprintf('============================================================\n\n');
fprintf('Filter RMSE Summary (Scenario S1 — Nominal):\n');
fprintf('  EKF : %.4f deg\n', EKF.rmse_deg);
fprintf('  UKF : %.4f deg\n', UKF.rmse_deg);
fprintf('  PF  : %.4f deg\n\n', PF.rmse_deg);
fprintf('Computation times:\n');
fprintf('  EKF : %.3f s\n', EKF.elapsed_s);
fprintf('  UKF : %.3f s\n', UKF.elapsed_s);
fprintf('  PF  : %.3f s  (%d particles)\n\n', PF.elapsed_s, PF.N_particles);
fprintf('Stress test RMSE table saved in STRESS.rmse_table\n');
fprintf('============================================================\n');


