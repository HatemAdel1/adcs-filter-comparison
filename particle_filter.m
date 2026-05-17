% =========================================================================
% particle_filter.m
% ADCS Filter Comparison Project — Module 2c: Particle Filter
%
% Sequential Importance Resampling (SIR) Particle Filter
% N_particles = 500  (balance accuracy vs computation)
%
% State : x = [q1; q2; q3; q4; bx; by; bz]  (7×1 per particle)
%
% Prediction : propagate each particle through dynamics + noise
% Update     : weight by Star Tracker likelihood (Gaussian)
% Resampling : systematic resampling when N_eff < N/2
%
% Inputs  : ENV, SENS structs
% Outputs : PF struct with estimated trajectory and error metrics
%
% Author  : Hatem Adel Saber
% =========================================================================

if ~exist('ENV','var') || ~exist('SENS','var')
    error('[PF] Run M0 and M1 first.');
end

fprintf('[PF] Running Particle Filter...\n');
tic;

N           = ENV.N;
dt          = ENV.dt;
N_part      = 500;     % number of particles
R_meas      = ENV.R_meas;

% -------------------------------------------------------------------------
% 1.  INITIALISE PARTICLES
% -------------------------------------------------------------------------
init_err = deg2rad(12);
init_ax  = [1;1;1]/norm([1;1;1]);
dq_init  = [init_ax*sin(init_err/2); cos(init_err/2)];
q0       = quat_mult(dq_init, ENV.q_true(:,1));

% Spread particles around initial estimate
particles = zeros(7, N_part);
for i = 1:N_part
    % Small random attitude perturbation
    dangle = 0.5*deg2rad(5) * randn(3,1);
    half   = norm(dangle)/2;
    if half < 1e-10
        dq_i = [dangle; 1];
    else
        dq_i = [dangle/norm(dangle)*sin(half); cos(half)];
    end
    q_i = quat_mult(dq_i, q0);
    b_i = 1e-4 * randn(3,1);    % small random initial bias
    particles(:,i) = [q_i; b_i];
end

weights = ones(1, N_part) / N_part;

% Storage
x_hist      = zeros(7, N);
att_err_deg = zeros(1, N);

% Initial estimate
x_hist(:,1) = weighted_mean(particles, weights);
x_hist(1:4,1) = x_hist(1:4,1) / norm(x_hist(1:4,1));

% -------------------------------------------------------------------------
% 2.  FILTER LOOP
% -------------------------------------------------------------------------
sigma_gyro = ENV.sigma_gyro * 2;   % process noise for particle spread
sigma_bias = ENV.sigma_bias * 2;

for k = 1:N-1

    % ---- PREDICT : propagate each particle ----------------------------
    for i = 1:N_part
        q_i = particles(1:4,i);
        b_i = particles(5:7,i);

        % Add process noise to gyro measurement
        omega_noisy = SENS.gyro_meas(:,k) - b_i + sigma_gyro*randn(3,1);

        % Integrate quaternion kinematics
        Xi    = [q_i(4)*eye(3)+skew(q_i(1:3)); -q_i(1:3)'];
        q_new = q_i + 0.5*dt*Xi*omega_noisy;
        q_new = q_new / norm(q_new);

        % Bias random walk
        b_new = b_i + sigma_bias*sqrt(dt)*randn(3,1);

        particles(:,i) = [q_new; b_new];
    end

    % ---- UPDATE : weight by likelihood --------------------------------
    if SENS.st_valid(k+1)
        z_meas = SENS.st_meas(:,k+1);   % 4×1 quaternion

        for i = 1:N_part
            q_i = particles(1:4,i);

            % Innovation as rotation vector (3×1)
            dq_i = quat_mult(z_meas, quat_inv(q_i));
            if dq_i(4) < 0, dq_i = -dq_i; end
            innov_i = 2*dq_i(1:3);

            % Gaussian likelihood  p(z|x_i)
            log_w = -0.5 * innov_i' * (R_meas \ innov_i);
            weights(i) = weights(i) * exp(log_w);
        end

        % Normalise weights
        w_sum = sum(weights);
        if w_sum < 1e-300 || isnan(w_sum)
            weights = ones(1,N_part) / N_part;   % weight collapse recovery
        else
            weights = weights / w_sum;
        end

        % Effective sample size
        N_eff = 1 / sum(weights.^2);

        % Systematic resampling when N_eff < N/2
        if N_eff < N_part/2
            particles = systematic_resample(particles, weights, N_part);
            weights   = ones(1,N_part) / N_part;
        end
    end

    % ---- ESTIMATE : weighted mean ------------------------------------
    x_est = weighted_mean(particles, weights);
    x_est(1:4) = x_est(1:4) / norm(x_est(1:4));
    x_hist(:,k+1) = x_est;
end

% -------------------------------------------------------------------------
% 3.  ATTITUDE ERROR
% -------------------------------------------------------------------------
for k = 1:N
    dq = quat_mult(quat_inv(ENV.q_true(:,k)), x_hist(1:4,k));
    if dq(4) < 0, dq = -dq; end
    att_err_deg(k) = 2*asin(min(1,norm(dq(1:3))))*(180/pi);
end

elapsed = toc;

% -------------------------------------------------------------------------
% 4.  PACKAGE RESULTS
% -------------------------------------------------------------------------
PF.x_hist       = x_hist;
PF.att_err_deg  = att_err_deg;
PF.rmse_deg     = sqrt(mean(att_err_deg(50:end).^2));
PF.elapsed_s    = elapsed;
PF.N_particles  = N_part;
PF.label        = 'PF';

fprintf('[PF] Done in %.3f s  (%d particles)\n', elapsed, N_part);
fprintf('[PF] RMSE (steady state): %.4f deg\n', PF.rmse_deg);

% -------------------------------------------------------------------------
% 5.  QUICK PLOT
% -------------------------------------------------------------------------
figure('Name','M2c — PF Results','NumberTitle','off','Position',[200 100 1100 420]);

subplot(1,2,1);
plot(ENV.t, att_err_deg, 'g', 'LineWidth', 1.0);
xlabel('Time [s]'); ylabel('Attitude error [deg]');
title('PF — Attitude estimation error');
yline(PF.rmse_deg,'k--',sprintf('RMSE = %.4f°',PF.rmse_deg));
grid on;

subplot(1,2,2);
bias_est = x_hist(5:7,:);
for ax = 1:3
    plot(ENV.t, rad2deg(bias_est(ax,:)-SENS.gyro_bias(ax,:))*3600); hold on;
end
xlabel('Time [s]'); ylabel('Bias error [arcsec/s]');
title('PF — Gyro bias estimation error');
legend('b_x','b_y','b_z'); grid on;

sgtitle('Module 2c — Particle Filter (500 particles)');
fprintf('[PF] PF struct ready.\n\n');

% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function x_mean = weighted_mean(particles, weights)
% Weighted mean — handles quaternion averaging via simple weighted sum + renorm
x_mean = sum(bsxfun(@times, particles, weights), 2);
end

function new_particles = systematic_resample(particles, weights, N)
% Systematic resampling — O(N), low variance
positions = ((0:N-1) + rand()) / N;
cumW      = cumsum(weights);
i = 1; j = 1;
new_particles = zeros(size(particles));
while i <= N
    if positions(i) <= cumW(j)
        new_particles(:,i) = particles(:,j);
        i = i + 1;
    else
        j = j + 1;
        if j > N, j = N; end
    end
end
end

function q_out = quat_mult(p, q)
pv = p(1:3); ps = p(4); qv = q(1:3); qs = q(4);
q_out = [ps*qv+qs*pv+cross(pv,qv); ps*qs-dot(pv,qv)];
q_out = q_out / norm(q_out);
end

function q_inv = quat_inv(q)
q_inv = [-q(1:3); q(4)];
end

function S = skew(v)
S = [0,-v(3),v(2); v(3),0,-v(1); -v(2),v(1),0];
end
