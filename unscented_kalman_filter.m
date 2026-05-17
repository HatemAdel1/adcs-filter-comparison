% =========================================================================
% unscented_kalman_filter.m
% ADCS Filter Comparison Project — Module 2b: Unscented Kalman Filter
%
% Uses sigma points — no Jacobians needed.
%
% State representation:
%   Nominal state : q(4×1) + bias(3×1)  — kept separately
%   Error state   : [d_phi(3×1); d_b(3×1)] = 6×1
%                   (avoids quaternion constraint singularity in P)
%
% Covariance P is always 6×6. Sigma points are generated in error space
% and composed onto the nominal quaternion before propagation.
%
% UKF parameters (Wan & Merwe):
%   alpha = 1e-3,  beta = 2,  kappa = 0
%
% Inputs  : ENV, SENS structs
% Outputs : UKF struct with estimated trajectory and error metrics
%
% Author  : Hatem Adel Saber
% MSc     : Space Electronics & Communications — Beni Suef University
% =========================================================================

if ~exist('ENV','var') || ~exist('SENS','var')
    error('[UKF] Run shared_simulation_environment and imu_and_star_tracker_simulation first.');
end

fprintf('[UKF] Running Unscented Kalman Filter (6-state error formulation)...\n');
tic;

N  = ENV.N;
dt = ENV.dt;
n  = 6;     % error-state dimension: 3 rotation + 3 bias

% Process noise (6×6 error-state)
sigma_gyro = ENV.sigma_gyro;
sigma_bias = ENV.sigma_bias;
q_phi = (sigma_gyro * dt)^2;
q_b   = (sigma_bias * dt)^2;
Q6    = diag([q_phi*ones(1,3), q_b*ones(1,3)]);

R_meas = ENV.R_meas;   % 3×3

% -------------------------------------------------------------------------
% 1.  UKF WEIGHTS
% -------------------------------------------------------------------------
alpha  = 1e-3;
beta   = 2;
kappa  = 0;
lambda = alpha^2 * (n + kappa) - n;
gamma  = sqrt(n + lambda);

Wm = [lambda/(n+lambda),  repmat(1/(2*(n+lambda)), 1, 2*n)];
Wc = [(lambda/(n+lambda) + 1 - alpha^2 + beta), repmat(1/(2*(n+lambda)), 1, 2*n)];

% -------------------------------------------------------------------------
% 2.  INITIALISATION
% -------------------------------------------------------------------------
init_err = deg2rad(12);
init_ax  = [1;1;1] / norm([1;1;1]);
dq_init  = [init_ax*sin(init_err/2); cos(init_err/2)];
q_est    = quat_mult(dq_init, ENV.q_true(:,1));
b_est    = zeros(3,1);

P = diag([1e-4*ones(1,3), 1e-6*ones(1,3)]);   % 6×6

x_hist = zeros(7, N);
P_hist = zeros(6, 6, N);
x_hist(:,1)   = [q_est; b_est];
P_hist(:,:,1) = P;

% -------------------------------------------------------------------------
% 3.  FILTER LOOP
% -------------------------------------------------------------------------
for k = 1:N-1

    % ============================================================
    % PREDICT
    % ============================================================

    % Sigma points in error space: 6×(2n+1)
    try
        Sc = gamma * chol(P, 'lower');
    catch
        P  = nearestSPD(P);
        Sc = gamma * chol(P, 'lower');
    end
    eps_sigma = [zeros(6,1), Sc, -Sc];   % 6×(2n+1)

    % Lift error sigma points to full state (q, b)
    q_sigma = zeros(4, 2*n+1);
    b_sigma = zeros(3, 2*n+1);
    for i = 1:2*n+1
        dq_i         = rotvec_to_quat(eps_sigma(1:3, i));
        q_sigma(:,i) = quat_mult(dq_i, q_est);
        b_sigma(:,i) = b_est + eps_sigma(4:6, i);
    end

    % Propagate each sigma point
    q_sp = zeros(4, 2*n+1);
    b_sp = zeros(3, 2*n+1);
    for i = 1:2*n+1
        [q_sp(:,i), b_sp(:,i)] = propagate(q_sigma(:,i), b_sigma(:,i), ...
                                             SENS.gyro_meas(:,k), dt);
    end

    % Predicted mean
    q_pred = sum(bsxfun(@times, q_sp, Wm), 2);
    q_pred = q_pred / norm(q_pred);
    b_pred = sum(bsxfun(@times, b_sp, Wm), 2);

    % Predicted covariance (6×6 in error space)
    P_pred = Q6;
    for i = 1:2*n+1
        dq_i   = quat_mult(q_sp(:,i), quat_inv(q_pred));
        if dq_i(4) < 0, dq_i = -dq_i; end
        eps_i  = [2*dq_i(1:3); b_sp(:,i) - b_pred];   % 6×1
        P_pred = P_pred + Wc(i) * (eps_i * eps_i');
    end
    P_pred = (P_pred + P_pred') / 2;

    % ============================================================
    % UPDATE
    % ============================================================
    if SENS.st_valid(k+1)

        try
            Sc2 = gamma * chol(P_pred, 'lower');
        catch
            P_pred = nearestSPD(P_pred);
            Sc2    = gamma * chol(P_pred, 'lower');
        end
        eps2 = [zeros(6,1), Sc2, -Sc2];

        q_su = zeros(4, 2*n+1);
        b_su = zeros(3, 2*n+1);
        for i = 1:2*n+1
            dq_i         = rotvec_to_quat(eps2(1:3,i));
            q_su(:,i)    = quat_mult(dq_i, q_pred);
            b_su(:,i)    = b_pred + eps2(4:6,i);
        end

        % Measurement sigma points (3×(2n+1))
        z_sigma = zeros(3, 2*n+1);
        for i = 1:2*n+1
            dq_i         = quat_mult(q_su(:,i), quat_inv(q_pred));
            if dq_i(4) < 0, dq_i = -dq_i; end
            z_sigma(:,i) = 2 * dq_i(1:3);
        end
        z_pred = sum(bsxfun(@times, z_sigma, Wm), 2);   % 3×1

        % Innovation and cross-covariances
        Pzz = R_meas;
        Pxz = zeros(6, 3);
        for i = 1:2*n+1
            dz_i = z_sigma(:,i) - z_pred;
            dq_i = quat_mult(q_su(:,i), quat_inv(q_pred));
            if dq_i(4) < 0, dq_i = -dq_i; end
            eps_i = [2*dq_i(1:3); b_su(:,i) - b_pred];   % 6×1
            Pzz   = Pzz + Wc(i) * (dz_i * dz_i');
            Pxz   = Pxz + Wc(i) * (eps_i * dz_i');
        end

        % Actual innovation
        z_meas  = SENS.st_meas(:,k+1);
        dq_meas = quat_mult(z_meas, quat_inv(q_pred));
        if dq_meas(4) < 0, dq_meas = -dq_meas; end
        innov   = 2*dq_meas(1:3) - z_pred;   % 3×1

        % Kalman gain and correction
        K      = Pxz / Pzz;      % 6×3
        dx     = K * innov;       % 6×1

        dq_cor = rotvec_to_quat(dx(1:3));
        q_est  = quat_mult(dq_cor, q_pred);
        b_est  = b_pred + dx(4:6);

        % Joseph form covariance update
        IKPxzT = eye(6) - K * (Pxz' / P_pred);
        P = IKPxzT * P_pred * IKPxzT' + K * R_meas * K';
    else
        q_est = q_pred;
        b_est = b_pred;
        P     = P_pred;
    end

    q_est = q_est / norm(q_est);
    P     = (P + P') / 2;

    x_hist(:,k+1)   = [q_est; b_est];
    P_hist(:,:,k+1) = P;
end

% -------------------------------------------------------------------------
% 4.  ATTITUDE ERROR & PACKAGING
% -------------------------------------------------------------------------
att_err_deg = zeros(1, N);
for k = 1:N
    dq = quat_mult(quat_inv(ENV.q_true(:,k)), x_hist(1:4,k));
    if dq(4) < 0, dq = -dq; end
    att_err_deg(k) = 2 * asin(min(1, norm(dq(1:3)))) * (180/pi);
end

elapsed = toc;

UKF.x_hist      = x_hist;
UKF.P_hist      = P_hist;
UKF.att_err_deg = att_err_deg;
UKF.rmse_deg    = sqrt(mean(att_err_deg(50:end).^2));
UKF.elapsed_s   = elapsed;
UKF.label       = 'UKF';

fprintf('[UKF] Done in %.3f s\n', elapsed);
fprintf('[UKF] RMSE (steady state): %.4f deg\n', UKF.rmse_deg);

% -------------------------------------------------------------------------
% 5.  QUICK PLOT
% -------------------------------------------------------------------------
figure('Name','UKF Results','NumberTitle','off','Position',[150 100 1100 420]);

subplot(1,2,1);
plot(ENV.t, att_err_deg, 'r', 'LineWidth', 1.0);
xlabel('Time [s]'); ylabel('Attitude error [deg]');
title('UKF — Attitude estimation error');
yline(UKF.rmse_deg,'k--',sprintf('RMSE = %.4f°', UKF.rmse_deg));
grid on;

subplot(1,2,2);
bias_est = x_hist(5:7,:);
for ax = 1:3
    plot(ENV.t, rad2deg(bias_est(ax,:) - SENS.gyro_bias(ax,:)) * 3600); hold on;
end
xlabel('Time [s]'); ylabel('Bias error [arcsec/s]');
title('UKF — Gyro bias estimation error');
legend('b_x','b_y','b_z'); grid on;

sgtitle('UKF — Unscented Kalman Filter (6-state error formulation)');
fprintf('[UKF] UKF struct ready for comparison.\n\n');

% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function [q_new, b_new] = propagate(q, b, omega_meas, dt)
omega_corr = omega_meas - b;
Xi    = [q(4)*eye(3) + skew(q(1:3)); -q(1:3)'];
q_new = q + 0.5 * dt * Xi * omega_corr;
q_new = q_new / norm(q_new);
b_new = b;
end

function dq = rotvec_to_quat(phi)
% Rotation vector phi (3×1) → unit quaternion [vec; scalar]
phi  = phi(:);
angle = norm(phi);
if angle < 1e-10
    dq = [phi * 0.5; 1];
else
    dq = [phi/angle * sin(angle/2); cos(angle/2)];
end
dq = dq / norm(dq);
end

function q_out = quat_mult(p, q)
pv = p(1:3); ps = p(4);
qv = q(1:3); qs = q(4);
q_out = [ps*qv(:) + qs*pv(:) + cross(pv(:), qv(:)); ps*qs - dot(pv,qv)];
q_out = q_out / norm(q_out);
end

function q_inv = quat_inv(q)
q_inv = [-q(1:3); q(4)];
end

function S = skew(v)
S = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end

function Aout = nearestSPD(A)
B    = (A + A') / 2;
[V,D] = eig(B);
d    = max(diag(D), 1e-10);
Aout = V * diag(d) * V';
Aout = (Aout + Aout') / 2;
end