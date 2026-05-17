% =========================================================================
% shared_simulation_environment.m
% ADCS Filter Comparison Project — Module 0: Shared Simulation Environment
%
% LEO CubeSat — Nadir Pointing Mission
% Altitude : 550 km circular orbit
% Sim time : 600 s  |  dt = 0.1 s
%
% Output   : struct ENV containing true attitude trajectory,
%             angular velocity, orbital parameters, and noise configs
%             for all three filters (EKF, UKF, PF).
%
% Author   : Hatem Adel Saber
% MSc      : Space Electronics & Communications — Beni Suef University
% =========================================================================

clear; clc; close all;
rng(42);   % fixed seed — identical noise for all three filters

% -------------------------------------------------------------------------
% 1.  SIMULATION TIME VECTOR
% -------------------------------------------------------------------------
dt      = 0.1;          % [s]  time step
t_end   = 600;          % [s]  total simulation time
t       = 0:dt:t_end;   % time vector
N       = length(t);    % number of steps

% -------------------------------------------------------------------------
% 2.  ORBITAL PARAMETERS  (LEO, 550 km circular)
% -------------------------------------------------------------------------
R_earth = 6371e3;       % [m]   Earth mean radius
mu      = 3.986e14;     % [m^3/s^2] Earth gravitational parameter
alt     = 550e3;        % [m]   orbital altitude
r_orb   = R_earth + alt;% [m]   orbital radius

n_orb   = sqrt(mu / r_orb^3);   % [rad/s] mean motion (orbital rate)
T_orb   = 2*pi / n_orb;         % [s]     orbital period (~5740 s)

fprintf('=== Orbital parameters ===\n');
fprintf('  Altitude     : %.0f km\n', alt/1e3);
fprintf('  Orbital rate : %.6f rad/s\n', n_orb);
fprintf('  Orbital period: %.1f s (%.2f min)\n', T_orb, T_orb/60);

% -------------------------------------------------------------------------
% 3.  TRUE ANGULAR VELOCITY  (nadir-pointing + small wobble)
%
%   omega_true = [n_orb*sin(small wobble); n_orb (pitch to track nadir);
%                  small roll oscillation]
%   All three axes have gentle sinusoidal variation representing realistic
%   attitude motion around the nadir-pointing equilibrium.
% -------------------------------------------------------------------------
omega_true = zeros(3, N);
for k = 1:N
    omega_true(1,k) =  0.5*n_orb * sin(2*pi*t(k)/120);   % roll
    omega_true(2,k) =  n_orb;                              % pitch (nadir)
    omega_true(3,k) =  0.3*n_orb * cos(2*pi*t(k)/180);   % yaw
end

% -------------------------------------------------------------------------
% 4.  TRUE ATTITUDE TRAJECTORY  (quaternion integration)
%
%   Quaternion convention : q = [q1; q2; q3; q4] with q4 = scalar part
%   Kinematics            : q_dot = 0.5 * Xi(q) * omega
% -------------------------------------------------------------------------

% initial attitude — 12 deg error (realistic post-detumbling)
angle0  = deg2rad(12);
axis0   = [1; 1; 1] / norm([1; 1; 1]);
q_true  = zeros(4, N);
q_true(:,1) = [axis0 * sin(angle0/2); cos(angle0/2)];   % [vector; scalar]

for k = 1:N-1
    q  = q_true(:,k);
    w  = omega_true(:,k);

    % Xi matrix  (4x3) : q_dot = 0.5 * Xi * omega
    Xi = [ q(4)*eye(3) + skew(q(1:3)); ...
          -q(1:3)' ];

    q_dot = 0.5 * Xi * w;

    % Euler integration + renormalize
    q_next = q + q_dot * dt;
    q_true(:,k+1) = q_next / norm(q_next);
end

% -------------------------------------------------------------------------
% 5.  NOISE PARAMETERS  (Nominal — Scenario S1)
% -------------------------------------------------------------------------

% --- IMU / Gyroscope ---
sigma_gyro   = 5e-5;     % [rad/s]   per-axis measurement noise std
sigma_bias   = 1e-6;     % [rad/s^2] bias drift rate std (random walk)

% --- Star Tracker ---
sigma_st_cb  = 48e-6;    % [rad]   cross-boresight  (10 arcsec)
sigma_st_bs  = 145e-6;   % [rad]   boresight / roll (30 arcsec)
f_st         = 1.0;      % [Hz]    update rate
step_st      = round(1 / (f_st * dt));   % update every this many steps

% --- Process noise covariance  Q  (7x7)  state = [q(4); bias(3)] ---
q_q   = (sigma_gyro * dt)^2;    % quaternion process noise variance
q_b   = (sigma_bias * dt)^2;    % bias process noise variance
Q     = diag([q_q*ones(1,4), q_b*ones(1,3)]);

% --- Measurement noise covariance  R  (3x3) ---
R_meas = diag([sigma_st_cb^2, sigma_st_cb^2, sigma_st_bs^2]);

% -------------------------------------------------------------------------
% 6.  PACKAGE INTO ENV STRUCT
% -------------------------------------------------------------------------
ENV.t           = t;
ENV.N           = N;
ENV.dt          = dt;

ENV.q_true      = q_true;
ENV.omega_true  = omega_true;

ENV.n_orb       = n_orb;
ENV.T_orb       = T_orb;
ENV.alt         = alt;

% noise — nominal (S1)
ENV.sigma_gyro  = sigma_gyro;
ENV.sigma_bias  = sigma_bias;
ENV.sigma_st_cb = sigma_st_cb;
ENV.sigma_st_bs = sigma_st_bs;
ENV.f_st        = f_st;
ENV.step_st     = step_st;
ENV.Q           = Q;
ENV.R_meas      = R_meas;

% stress-test multipliers (used by Module 5)
ENV.stress.high_imu_mult  = 10;    % ×10 gyro noise
ENV.stress.high_st_mult   = 6;     % ×6  star tracker noise
ENV.stress.outage_dur     = 100;   % steps blind per 600 steps
ENV.stress.outage_period  = 600;   % steps

fprintf('\n=== Simulation parameters ===\n');
fprintf('  Steps        : %d\n', N);
fprintf('  dt           : %.2f s\n', dt);
fprintf('  sigma_gyro   : %.2e rad/s\n', sigma_gyro);
fprintf('  sigma_st_cb  : %.2e rad (%.1f arcsec)\n', sigma_st_cb, sigma_st_cb*206265);
fprintf('  ST update    : every %d steps (%.1f Hz)\n', step_st, f_st);
fprintf('  Q size       : %dx%d\n', size(Q));
fprintf('  R_meas size  : %dx%d\n', size(R_meas));

% -------------------------------------------------------------------------
% 7.  QUICK SANITY PLOT — True Attitude & Angular Velocity
% -------------------------------------------------------------------------
figure('Name','M0 — True Attitude & Angular Velocity','NumberTitle','off', ...
       'Position',[100 100 1100 500]);

subplot(1,2,1);
plot(t, q_true(1,:), 'b',  t, q_true(2,:), 'r', ...
     t, q_true(3,:), 'g',  t, q_true(4,:), 'k', 'LineWidth',1.2);
xlabel('Time [s]'); ylabel('Quaternion components');
title('True attitude — q(t)');
legend('q_1','q_2','q_3','q_4','Location','best');
grid on;

subplot(1,2,2);
plot(t, rad2deg(omega_true(1,:)), 'b', ...
     t, rad2deg(omega_true(2,:)), 'r', ...
     t, rad2deg(omega_true(3,:)), 'g', 'LineWidth',1.2);
xlabel('Time [s]'); ylabel('\omega [deg/s]');
title('True angular velocity — \omega(t)');
legend('\omega_x','\omega_y','\omega_z','Location','best');
grid on;

sgtitle('Module 0 — Shared Environment (LEO CubeSat, Nadir Pointing)');

fprintf('\n[M0] Done. ENV struct ready with %d fields.\n', length(fieldnames(ENV)));
fprintf('[M0] Save ENV before running M1:\n');
fprintf('       save(''ENV.mat'', ''ENV'')\n\n');

% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function S = skew(v)
% skew(v)  —  3x3 skew-symmetric matrix of vector v
S = [  0    -v(3)  v(2); ...
       v(3)  0    -v(1); ...
      -v(2)  v(1)  0   ];
end
