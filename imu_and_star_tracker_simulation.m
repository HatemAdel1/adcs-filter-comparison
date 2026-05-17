% =========================================================================
% imu_and_star_tracker_simulation.m
% ADCS Filter Comparison Project — Module 1: Sensor Simulation
%
% Generates noisy IMU (gyroscope) measurements and Star Tracker
% attitude observations from the true trajectory in ENV.
%
% Inputs  : ENV struct from shared_simulation_environment.m
% Outputs : SENS struct with gyro_meas, st_meas, st_valid flags
%
% Author  : Hatem Adel Saber
% =========================================================================

% Load environment (run M0 first, then save ENV)
% load('ENV.mat');      % uncomment when running standalone
% If running after M0 in same session, ENV is already in workspace.

if ~exist('ENV','var')
    error('[M1] ENV struct not found. Run shared_simulation_environment.m first.');
end

fprintf('[M1] Generating sensor measurements...\n');

N          = ENV.N;
dt         = ENV.dt;
q_true     = ENV.q_true;
omega_true = ENV.omega_true;

% -------------------------------------------------------------------------
% 1.  IMU — GYROSCOPE SIMULATION
%
%   Model : omega_meas = omega_true + bias + white_noise
%           bias(k+1)  = bias(k)   + sigma_bias * sqrt(dt) * randn
% -------------------------------------------------------------------------
sigma_gyro = ENV.sigma_gyro;
sigma_bias = ENV.sigma_bias;

gyro_bias  = zeros(3, N);
gyro_meas  = zeros(3, N);

% initial bias — small random offset (realistic)
gyro_bias(:,1) = sigma_gyro * 2 * randn(3,1);

for k = 1:N-1
    % bias random walk
    gyro_bias(:,k+1) = gyro_bias(:,k) + sigma_bias * sqrt(dt) * randn(3,1);

    % measurement = truth + bias + white noise
    gyro_meas(:,k) = omega_true(:,k) + gyro_bias(:,k) ...
                     + sigma_gyro * randn(3,1);
end
gyro_meas(:,N) = omega_true(:,N) + gyro_bias(:,N) + sigma_gyro * randn(3,1);

% -------------------------------------------------------------------------
% 2.  STAR TRACKER — ATTITUDE OBSERVATION SIMULATION
%
%   Outputs a noisy quaternion measurement at 1 Hz (every step_st steps).
%   Noise model : small-angle rotation error added to true quaternion,
%                 then renormalized.
%   Error vector: [err_cb; err_cb; err_bs]  (cross-boresight x2, boresight)
% -------------------------------------------------------------------------
sigma_cb  = ENV.sigma_st_cb;
sigma_bs  = ENV.sigma_st_bs;
step_st   = ENV.step_st;

st_meas   = zeros(4, N);   % quaternion measurement (only valid steps used)
st_valid  = false(1, N);   % true when ST update is available

for k = 1:N
    if mod(k-1, step_st) == 0    % ST fires at these steps
        % small-angle attitude error in body frame
        d_angle = [sigma_cb * randn; ...
                   sigma_cb * randn; ...
                   sigma_bs * randn];

        half_angle = norm(d_angle) / 2;
        if half_angle < 1e-10
            dq = [0; 0; 0; 1];
        else
            dq = [d_angle/norm(d_angle) * sin(half_angle); cos(half_angle)];
        end

        % compose: q_meas = dq ⊗ q_true
        st_meas(:,k) = quat_mult(dq, q_true(:,k));
        st_valid(k)  = true;
    end
end

% -------------------------------------------------------------------------
% 3.  PACKAGE INTO SENS STRUCT
% -------------------------------------------------------------------------
SENS.gyro_meas  = gyro_meas;
SENS.gyro_bias  = gyro_bias;   % ground truth bias (for error analysis)
SENS.st_meas    = st_meas;
SENS.st_valid   = st_valid;
SENS.N_st_meas  = sum(st_valid);

fprintf('[M1] Gyroscope : %d samples  |  sigma = %.2e rad/s\n', N, sigma_gyro);
fprintf('[M1] Star Tracker: %d updates  |  sigma_cb = %.1f arcsec\n', ...
        sum(st_valid), sigma_cb*206265);

% -------------------------------------------------------------------------
% 4.  QUICK SANITY PLOTS
% -------------------------------------------------------------------------
t = ENV.t;

figure('Name','M1 — Sensor Measurements','NumberTitle','off', ...
       'Position',[100 100 1200 550]);

subplot(2,2,1);
plot(t, rad2deg(gyro_meas(1,:)), 'b', 'LineWidth', 0.8); hold on;
plot(t, rad2deg(omega_true(1,:)), 'k--', 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('[deg/s]');
title('Gyro x — measured vs true');
legend('Measured','True'); grid on;

subplot(2,2,2);
bias_err = gyro_bias - 0;   % true bias trajectory
plot(t, rad2deg(gyro_bias(1,:)), 'r', ...
     t, rad2deg(gyro_bias(2,:)), 'g', ...
     t, rad2deg(gyro_bias(3,:)), 'b', 'LineWidth', 0.9);
xlabel('Time [s]'); ylabel('[deg/s]');
title('True gyro bias — random walk');
legend('b_x','b_y','b_z'); grid on;

subplot(2,2,3);
st_t = t(st_valid);
err_q = zeros(1, sum(st_valid));
idx = 0;
for k = find(st_valid)
    idx = idx + 1;
    dq_err = quat_mult(quat_inv(q_true(:,k)), st_meas(:,k));
    err_q(idx) = 2 * asin(norm(dq_err(1:3))) * (180/pi) * 3600; % arcsec
end
plot(st_t, err_q, 'mo', 'MarkerSize', 3, 'MarkerFaceColor','m');
xlabel('Time [s]'); ylabel('Error [arcsec]');
title('Star Tracker attitude error');
yline(10*ones(1,1),'k--','10 arcsec spec'); grid on;

subplot(2,2,4);
gyro_err_x = rad2deg((gyro_meas(1,:) - omega_true(1,:))) * 3600;
histogram(gyro_err_x, 40, 'FaceColor','b', 'EdgeColor','none', 'FaceAlpha',0.7);
xlabel('Error [arcsec/s equiv]'); ylabel('Count');
title('Gyro noise distribution (x-axis)');
grid on;

sgtitle('Module 1 — IMU & Star Tracker Simulation');

fprintf('[M1] Done. SENS struct ready.\n');
fprintf('[M1] Save: save(''SENS.mat'',''SENS'')  or keep in workspace for M2.\n\n');

% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function q_out = quat_mult(p, q)
% Quaternion multiplication: q_out = p ⊗ q
% Convention: q = [vector(3); scalar]
pv = p(1:3); ps = p(4);
qv = q(1:3); qs = q(4);
q_out = [ps*qv + qs*pv + cross(pv,qv); ...
         ps*qs - dot(pv,qv)];
q_out = q_out / norm(q_out);
end

function q_inv = quat_inv(q)
% Quaternion inverse (conjugate for unit quaternion)
q_inv = [-q(1:3); q(4)];
end
