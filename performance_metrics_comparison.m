% =========================================================================
% performance_metrics_comparison.m
% ADCS Filter Comparison Project — Module 3: Performance Metrics
%
% Compares EKF, UKF, and PF on:
%   1. RMSE (steady-state attitude error)
%   2. Convergence time
%   3. Computation time
%   4. Max peak error
%   5. 3-sigma bounds consistency
%
% Inputs  : ENV, EKF, UKF, PF structs
% Outputs : RESULTS struct + printed comparison table
%
% Author  : Hatem Adel Saber
% =========================================================================

if ~exist('EKF','var') || ~exist('UKF','var') || ~exist('PF','var')
    error('[M3] Run M2a, M2b, M2c first to generate EKF, UKF, PF structs.');
end

fprintf('\n========== MODULE 3: FILTER COMPARISON RESULTS ==========\n\n');

filters = {EKF, UKF, PF};
names   = {'EKF', 'UKF', 'PF'};
colors  = {'b',   'r',   'g'};
N       = ENV.N;
t       = ENV.t;

% -------------------------------------------------------------------------
% 1.  CONVERGENCE TIME  (first time error < 0.5 deg and stays < 0.5 deg)
% -------------------------------------------------------------------------
thresh = 0.5;   % [deg]
conv_time = zeros(1,3);
for f = 1:3
    err = filters{f}.att_err_deg;
    for k = 1:N-50
        if all(err(k:k+50) < thresh)
            conv_time(f) = t(k);
            break;
        end
    end
    if conv_time(f) == 0
        conv_time(f) = NaN;
    end
end

% -------------------------------------------------------------------------
% 2.  COMPILE METRICS TABLE
% -------------------------------------------------------------------------
RESULTS.filter_names  = names;
RESULTS.rmse_deg      = cellfun(@(f) f.rmse_deg,    filters);
RESULTS.peak_err_deg  = cellfun(@(f) max(f.att_err_deg(50:end)), filters);
RESULTS.conv_time_s   = conv_time;
RESULTS.elapsed_s     = cellfun(@(f) f.elapsed_s,   filters);

% Print formatted table
fprintf('%-20s  %-10s  %-10s  %-10s\n', 'Metric', 'EKF', 'UKF', 'PF');
fprintf('%s\n', repmat('-',1,56));
fprintf('%-20s  %-10.4f  %-10.4f  %-10.4f  [deg]\n', ...
        'RMSE (steady)', RESULTS.rmse_deg(1), RESULTS.rmse_deg(2), RESULTS.rmse_deg(3));
fprintf('%-20s  %-10.4f  %-10.4f  %-10.4f  [deg]\n', ...
        'Peak error', RESULTS.peak_err_deg(1), RESULTS.peak_err_deg(2), RESULTS.peak_err_deg(3));
fprintf('%-20s  %-10.1f  %-10.1f  %-10.1f  [s]\n', ...
        'Convergence time', RESULTS.conv_time_s(1), RESULTS.conv_time_s(2), RESULTS.conv_time_s(3));
fprintf('%-20s  %-10.3f  %-10.3f  %-10.3f  [s]\n', ...
        'Compute time', RESULTS.elapsed_s(1), RESULTS.elapsed_s(2), RESULTS.elapsed_s(3));

% Relative improvement UKF vs EKF
dRMSE_ukf = (RESULTS.rmse_deg(1) - RESULTS.rmse_deg(2)) / RESULTS.rmse_deg(1) * 100;
dRMSE_pf  = (RESULTS.rmse_deg(1) - RESULTS.rmse_deg(3)) / RESULTS.rmse_deg(1) * 100;
fprintf('\n  UKF RMSE improvement over EKF : %+.1f%%\n', dRMSE_ukf);
fprintf('  PF  RMSE improvement over EKF : %+.1f%%\n\n', dRMSE_pf);

% -------------------------------------------------------------------------
% 3.  MAIN COMPARISON PLOT  (all 3 filters on one figure)
% -------------------------------------------------------------------------
figure('Name','M3 — Filter Comparison','NumberTitle','off', ...
       'Position',[50 50 1400 700]);

% --- Subplot 1: Error trajectories ---
subplot(2,2,1);
hold on;
for f = 1:3
    plot(t, filters{f}.att_err_deg, colors{f}, 'LineWidth', 1.2);
end
xlabel('Time [s]'); ylabel('Attitude error [deg]');
title('Attitude estimation error — all filters');
legend('EKF','UKF','PF','Location','northeast');
yline(thresh,'k--','0.5° threshold');
grid on; hold off;

% --- Subplot 2: RMSE bar chart ---
subplot(2,2,2);
b = bar(RESULTS.rmse_deg, 0.5, 'FaceColor','flat');
b.CData = [0 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];
set(gca,'XTickLabel', names);
xlabel('Filter'); ylabel('RMSE [deg]');
title('Steady-state RMSE comparison');
text(1:3, RESULTS.rmse_deg + 0.0002, ...
     arrayfun(@(x)sprintf('%.4f°',x), RESULTS.rmse_deg,'UniformOutput',false), ...
     'HorizontalAlignment','center','FontSize',9);
grid on;

% --- Subplot 3: Computation time ---
subplot(2,2,3);
b2 = bar(RESULTS.elapsed_s, 0.5, 'FaceColor','flat');
b2.CData = [0 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];
set(gca,'XTickLabel', names);
xlabel('Filter'); ylabel('Time [s]');
title('Computation time');
text(1:3, RESULTS.elapsed_s + 0.01, ...
     arrayfun(@(x)sprintf('%.3fs',x), RESULTS.elapsed_s,'UniformOutput',false), ...
     'HorizontalAlignment','center','FontSize',9);
grid on;

% --- Subplot 4: Error histogram ---
subplot(2,2,4);
hold on;
for f = 1:3
    err_ss = filters{f}.att_err_deg(50:end);
    histogram(err_ss, 40, 'FaceAlpha', 0.5, 'EdgeColor','none', ...
              'DisplayName', names{f});
end
xlabel('Attitude error [deg]'); ylabel('Count');
title('Error distribution (steady state)');
legend; grid on; hold off;

sgtitle('Module 3 — EKF vs UKF vs PF: Performance Comparison (S1 Nominal)');

fprintf('[M3] Comparison complete. RESULTS struct ready.\n');
fprintf('[M3] Proceed to M4 for advanced visualization.\n\n');
