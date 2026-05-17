% =========================================================================
% extended_kalman_filter.m
% ADCS Filter Comparison Project — Module 2a: Extended Kalman Filter
%
% State vector : x = [q1; q2; q3; q4; bx; by; bz]  (7×1)
%   q = attitude quaternion  [vector; scalar]
%   b = gyroscope bias       [rad/s]
%
% Prediction  : quaternion kinematics + bias random walk
% Update      : Star Tracker quaternion measurement
% Jacobians   : computed analytically
%
% Inputs  : ENV, SENS structs
% Outputs : EKF struct with estimated trajectory and error metrics
%
% Author  : Hatem Adel Saber
% =========================================================================

if ~exist('ENV','var') || ~exist('SENS','var')
    error('[EKF] Run M0 and M1 first to generate ENV and SENS.');
end

fprintf('[EKF] Running Extended Kalman Filter...\n');
tic;

N      = ENV.N;
dt     = ENV.dt;
Q      = ENV.Q;
R_meas = ENV.R_meas;

% -------------------------------------------------------------------------
% 1.  INITIALISATION
% -------------------------------------------------------------------------
% Initial state — perturb true quaternion slightly (filter doesn't know truth)
init_err_angle = deg2rad(12);
init_axis      = [1;1;1]/norm([1;1;1]);
dq_init = [init_axis*sin(init_err_angle/2); cos(init_err_angle/2)];
q0      = quat_mult(dq_init, ENV.q_true(:,1));

x_est  = [q0; zeros(3,1)];   % [q; bias=0]

% Initial covariance P
P = diag([1e-4*ones(1,4), 1e-6*ones(1,3)]);

% Storage
x_hist = zeros(7, N);
P_hist = zeros(7, 7, N);
x_hist(:,1)   = x_est;
P_hist(:,:,1) = P;

% -------------------------------------------------------------------------
% 2.  FILTER LOOP
% -------------------------------------------------------------------------
for k = 1:N-1

    % --- Extract current estimate ---
    q  = x_est(1:4);
    b  = x_est(5:7);

    % --- Corrected angular velocity ---
    omega_corr = SENS.gyro_meas(:,k) - b;

    % ---- PREDICT -------------------------------------------------------

    % State transition: integrate quaternion kinematics
    Xi    = [q(4)*eye(3) + skew(q(1:3)); -q(1:3)'];   % 4×3
    q_dot = 0.5 * Xi * omega_corr;
    q_new = q + q_dot * dt;
    q_new = q_new / norm(q_new);
    b_new = b;                   % bias modelled as constant between updates

    x_pred = [q_new; b_new];

    % Jacobian F = d(f)/dx  (7×7)
    F = ekf_jacobian_F(q, omega_corr, dt);

    % Predicted covariance
    P_pred = F * P * F' + Q;

    % ---- UPDATE (only when Star Tracker fires) --------------------------
    if SENS.st_valid(k+1)
        q_p   = x_pred(1:4);
        z_meas = SENS.st_meas(:,k+1);   % 4×1 quaternion measurement

        % Innovation (quaternion difference as rotation vector)
        dq_innov = quat_mult(z_meas, quat_inv(q_p));
        if dq_innov(4) < 0, dq_innov = -dq_innov; end
        innov = 2 * dq_innov(1:3);       % 3×1 small-angle approximation

        % Measurement Jacobian H (3×7) — maps state to measurement space
        H = ekf_jacobian_H(q_p);

        % Kalman gain
        S_innov = H * P_pred * H' + R_meas;   % 3×3
        K       = P_pred * H' / S_innov;       % 7×3

        % Corrected state
        dx      = K * innov;
        q_upd   = quat_plus_rotvec(x_pred(1:4), dx(1:4));
        b_upd   = x_pred(5:7) + dx(5:7);
        x_est   = [q_upd; b_upd];

        % Joseph form covariance update (numerically stable)
        IKH     = eye(7) - K * H;
        P       = IKH * P_pred * IKH' + K * R_meas * K';
    else
        x_est = x_pred;
        P     = P_pred;
    end

    % Enforce quaternion normalisation
    x_est(1:4) = x_est(1:4) / norm(x_est(1:4));

    x_hist(:,k+1)   = x_est;
    P_hist(:,:,k+1) = P;
end

% -------------------------------------------------------------------------
% 3.  COMPUTE ATTITUDE ERROR
% -------------------------------------------------------------------------
att_err_deg = zeros(1, N);
for k = 1:N
    q_est  = x_hist(1:4, k);
    q_true = ENV.q_true(:, k);
    dq     = quat_mult(quat_inv(q_true), q_est);
    if dq(4) < 0, dq = -dq; end
    att_err_deg(k) = 2 * asin(min(1, norm(dq(1:3)))) * (180/pi);
end

elapsed = toc;

% -------------------------------------------------------------------------
% 4.  PACKAGE RESULTS
% -------------------------------------------------------------------------
EKF.x_hist      = x_hist;
EKF.P_hist      = P_hist;
EKF.att_err_deg = att_err_deg;
EKF.rmse_deg    = sqrt(mean(att_err_deg(50:end).^2));   % skip transient
EKF.elapsed_s   = elapsed;
EKF.label       = 'EKF';

fprintf('[EKF] Done in %.3f s\n', elapsed);
fprintf('[EKF] RMSE (steady state): %.4f deg\n', EKF.rmse_deg);

% -------------------------------------------------------------------------
% 5.  QUICK PLOT
% -------------------------------------------------------------------------
figure('Name','M2a — EKF Results','NumberTitle','off','Position',[100 100 1100 420]);

subplot(1,2,1);
plot(ENV.t, att_err_deg, 'b', 'LineWidth', 1.0);
xlabel('Time [s]'); ylabel('Attitude error [deg]');
title('EKF — Attitude estimation error');
grid on;
yline(EKF.rmse_deg,'r--',sprintf('RMSE = %.4f°',EKF.rmse_deg));

subplot(1,2,2);
bias_est  = x_hist(5:7,:);
bias_true = SENS.gyro_bias;
for ax = 1:3
    plot(ENV.t, rad2deg(bias_est(ax,:) - bias_true(ax,:)) * 3600); hold on;
end
xlabel('Time [s]'); ylabel('Bias error [arcsec/s]');
title('EKF — Gyro bias estimation error');
legend('b_x','b_y','b_z'); grid on;

sgtitle('Module 2a — EKF');
fprintf('[EKF] EKF struct ready for Module 3 comparison.\n\n');

% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function F = ekf_jacobian_F(q, omega, dt)
% Linearised state transition Jacobian  F = I + df/dx * dt
% State: [q(4); bias(3)]
qv = q(1:3); qs = q(4);
W  = [-skew(omega) + 0*eye(3), -0.5*[qs*eye(3)+skew(qv)];
       omega',                   0,   0,   0             ];
% Build full 7×7 Jacobian analytically
Fqq = eye(4) + 0.5*dt*Omega_mat(omega);   % 4×4
Fqb = -0.5*dt*[qs*eye(3)+skew(qv); -qv']; % 4×3
Fbq = zeros(3,4);                          % 3×4
Fbb = eye(3);                              % 3×3
F   = [Fqq, Fqb; Fbq, Fbb];
end

function Om = Omega_mat(w)
% 4×4 Omega matrix for quaternion kinematics: q_dot = 0.5*Omega*q
Om = [  0,   w(3), -w(2),  w(1); ...
       -w(3), 0,    w(1),  w(2); ...
        w(2),-w(1), 0,     w(3); ...
       -w(1),-w(2),-w(3),  0   ];
end

function H = ekf_jacobian_H(q)
% Measurement Jacobian H (3×7)
% Innovation = 2*vec(dq) ≈ H * dx  around current estimate
% Linearisation of quaternion difference w.r.t. quaternion state
qs = q(4); qv = q(1:3);
% dinnov/dq  (3×4)
Hq = 2*[qs*eye(3) + skew(qv), -qv];  % simplified first-order
H  = [Hq, zeros(3,3)];               % bias doesn't affect ST measurement
end

function q_out = quat_plus_rotvec(q, dq_state)
% Add a small quaternion correction dq_state(1:4) to q
% dq_state(1:3) = small rotation vector, dq_state(4) = scalar correction
dq_vec = dq_state(1:3);
half   = norm(dq_vec)/2;
if half < 1e-10
    dq = [dq_vec; 1];
else
    dq = [dq_vec/norm(dq_vec)*sin(half); cos(half)];
end
q_out = quat_mult(dq, q);
q_out = q_out / norm(q_out);
end

function q_out = quat_mult(p, q)
pv = p(1:3); ps = p(4);
qv = q(1:3); qs = q(4);
q_out = [ps*qv + qs*pv + cross(pv,qv); ps*qs - dot(pv,qv)];
q_out = q_out / norm(q_out);
end

function q_inv = quat_inv(q)
q_inv = [-q(1:3); q(4)];
end

function S = skew(v)
S = [0,-v(3),v(2); v(3),0,-v(1); -v(2),v(1),0];
end
