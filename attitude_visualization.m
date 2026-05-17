% =========================================================================
% attitude_visualization.m
% ADCS Filter Comparison Project — Module 4: Advanced Visualization
%
% Produces publication-quality figures:
%   Fig 1 — Euler angle estimation (roll, pitch, yaw)
%   Fig 2 — 3D attitude trajectory (quaternion → Euler)
%   Fig 3 — Bias estimation comparison
%   Fig 4 — Filter uncertainty (3-sigma envelope from P matrix)
%
% Author  : Hatem Adel Saber
% =========================================================================

if ~exist('EKF','var') || ~exist('UKF','var') || ~exist('PF','var')
    error('[M4] Run M2a, M2b, M2c first.');
end

t  = ENV.t;
N  = ENV.N;

% -------------------------------------------------------------------------
% 1.  CONVERT QUATERNIONS TO EULER ANGLES  (ZYX convention)
% -------------------------------------------------------------------------
euler_true = zeros(3,N);
euler_ekf  = zeros(3,N);
euler_ukf  = zeros(3,N);
euler_pf   = zeros(3,N);

for k = 1:N
    euler_true(:,k) = quat2euler(ENV.q_true(:,k));
    euler_ekf(:,k)  = quat2euler(EKF.x_hist(1:4,k));
    euler_ukf(:,k)  = quat2euler(UKF.x_hist(1:4,k));
    euler_pf(:,k)   = quat2euler(PF.x_hist(1:4,k));
end

euler_true = rad2deg(euler_true);
euler_ekf  = rad2deg(euler_ekf);
euler_ukf  = rad2deg(euler_ukf);
euler_pf   = rad2deg(euler_pf);

% -------------------------------------------------------------------------
% 2.  FIG 1 — Euler angle estimation (3×1 subplots)
% -------------------------------------------------------------------------
angle_labels = {'Roll \phi [deg]', 'Pitch \theta [deg]', 'Yaw \psi [deg]'};

figure('Name','M4 — Euler Angle Estimation','NumberTitle','off', ...
       'Position',[50 50 1300 650]);

for ax = 1:3
    subplot(3,1,ax);
    plot(t, euler_true(ax,:), 'k', 'LineWidth', 1.5, 'DisplayName','True'); hold on;
    plot(t, euler_ekf(ax,:),  'b--','LineWidth', 1.0, 'DisplayName','EKF');
    plot(t, euler_ukf(ax,:),  'r:', 'LineWidth', 1.0, 'DisplayName','UKF');
    plot(t, euler_pf(ax,:),   'g-.','LineWidth', 1.0, 'DisplayName','PF');
    ylabel(angle_labels{ax});
    if ax == 1, title('Euler angle estimation — true vs estimated'); end
    if ax == 3, xlabel('Time [s]'); end
    legend('Location','best'); grid on; hold off;
end

% -------------------------------------------------------------------------
% 3.  FIG 2 — 3D attitude trajectory
% -------------------------------------------------------------------------
figure('Name','M4 — 3D Attitude Trajectory','NumberTitle','off', ...
       'Position',[100 100 900 700]);

plot3(euler_true(1,:), euler_true(2,:), euler_true(3,:), ...
      'k', 'LineWidth', 1.8, 'DisplayName','True');
hold on;
plot3(euler_ekf(1,:), euler_ekf(2,:), euler_ekf(3,:), ...
      'b--','LineWidth',1.0,'DisplayName','EKF');
plot3(euler_ukf(1,:), euler_ukf(2,:), euler_ukf(3,:), ...
      'r:', 'LineWidth',1.0,'DisplayName','UKF');
plot3(euler_pf(1,:),  euler_pf(2,:),  euler_pf(3,:), ...
      'g-.','LineWidth',1.0,'DisplayName','PF');

% Mark start and end
plot3(euler_true(1,1),  euler_true(2,1),  euler_true(3,1),  ...
      'ko','MarkerSize',8,'MarkerFaceColor','k');
plot3(euler_true(1,end),euler_true(2,end),euler_true(3,end), ...
      'ks','MarkerSize',8,'MarkerFaceColor','y');

xlabel('Roll \phi [deg]'); ylabel('Pitch \theta [deg]'); zlabel('Yaw \psi [deg]');
title('3D Attitude Trajectory — True vs Estimated (Roll-Pitch-Yaw)');
legend('Location','best'); grid on; view(35,25); hold off;

% -------------------------------------------------------------------------
% 4.  FIG 3 — Bias estimation comparison
% -------------------------------------------------------------------------
figure('Name','M4 — Bias Estimation','NumberTitle','off', ...
       'Position',[150 80 1300 500]);

bias_labels = {'b_x','b_y','b_z'};
for ax = 1:3
    subplot(1,3,ax);
    true_bias_deg = rad2deg(SENS.gyro_bias(ax,:)) * 3600;
    ekf_bias_deg  = rad2deg(EKF.x_hist(4+ax,:))  * 3600;
    ukf_bias_deg  = rad2deg(UKF.x_hist(4+ax,:))  * 3600;
    pf_bias_deg   = rad2deg(PF.x_hist(4+ax,:))   * 3600;

    plot(t, true_bias_deg,'k','LineWidth',1.5,'DisplayName','True'); hold on;
    plot(t, ekf_bias_deg, 'b--','LineWidth',0.9,'DisplayName','EKF');
    plot(t, ukf_bias_deg, 'r:', 'LineWidth',0.9,'DisplayName','UKF');
    plot(t, pf_bias_deg,  'g-.','LineWidth',0.9,'DisplayName','PF');
    xlabel('Time [s]'); ylabel('[arcsec/s]');
    title(sprintf('Gyro bias %s', bias_labels{ax}));
    legend('Location','best'); grid on; hold off;
end
sgtitle('Module 4 — Gyro Bias Estimation Comparison');

% -------------------------------------------------------------------------
% 5.  FIG 4 — EKF 3-sigma uncertainty envelope
% -------------------------------------------------------------------------
figure('Name','M4 — EKF Uncertainty Envelope','NumberTitle','off', ...
       'Position',[200 60 1000 400]);

% Extract 3-sigma from P diagonal (quaternion x-component as proxy)
sigma3_q1 = zeros(1,N);
for k = 1:N
    sigma3_q1(k) = 3 * sqrt(max(0, EKF.P_hist(1,1,k)));
end

err_q1 = EKF.x_hist(1,:) - ENV.q_true(1,:);
plot(t, err_q1,      'b',  'LineWidth', 1.0, 'DisplayName','q_1 error'); hold on;
plot(t,  sigma3_q1,  'r--','LineWidth', 0.8, 'DisplayName','+3\sigma');
plot(t, -sigma3_q1,  'r--','LineWidth', 0.8, 'DisplayName','-3\sigma');
fill([t, fliplr(t)], [sigma3_q1, fliplr(-sigma3_q1)], 'r', ...
     'FaceAlpha', 0.08, 'EdgeColor','none');

xlabel('Time [s]'); ylabel('q_1 error');
title('EKF — 3\sigma uncertainty envelope (q_1 component)');
legend('Location','best'); grid on; hold off;

fprintf('[M4] All visualization figures generated.\n');
fprintf('[M4] Ready for Module 5 stress tests.\n\n');

% =========================================================================
% LOCAL HELPER FUNCTION
% =========================================================================
function euler = quat2euler(q)
% ZYX convention: q = [qx; qy; qz; qw]
qx=q(1); qy=q(2); qz=q(3); qw=q(4);
roll  = atan2(2*(qw*qx+qy*qz), 1-2*(qx^2+qy^2));
pitch = asin(max(-1,min(1, 2*(qw*qy-qz*qx))));
yaw   = atan2(2*(qw*qz+qx*qy), 1-2*(qy^2+qz^2));
euler = [roll; pitch; yaw];
end
