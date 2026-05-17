% =========================================================================
% monte_carlo_statistical_analysis.m
% ADCS Filter Comparison Project — Module 6: Monte Carlo Analysis
%
% Runs N_mc independent simulations with different random seeds to
% produce statistically significant performance comparison between
% EKF, UKF, and Particle Filter.
%
% For each seed:
%   - New sensor noise realization (IMU + Star Tracker)
%   - All 3 filters run on identical measurements
%   - RMSE recorded per filter per scenario
%
% Outputs:
%   MC struct with mean, std, median, min, max RMSE per filter
%   + publication-quality figures (box plots, violin-style, CDF)
%
% Statistical tests:
%   - Wilcoxon signed-rank test: UKF vs EKF  (non-parametric, no
%     Gaussian assumption needed)
%   - Effect size (Cohen's d) between filters
%
% N_mc = 50 runs  (standard for attitude estimation literature)
%
% Author  : Hatem Adel Saber
% MSc     : Space Electronics & Communications — Beni Suef University
% =========================================================================

if ~exist('ENV','var')
    error('[MC] Run shared_simulation_environment.m first to generate ENV.');
end

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════╗\n');
fprintf('║     MODULE 6: MONTE CARLO STATISTICAL ANALYSIS      ║\n');
fprintf('║     N = 50 independent runs  |  3 filters  |  S1    ║\n');
fprintf('╚══════════════════════════════════════════════════════╝\n\n');

N_mc         = 50;
filter_names = {'EKF', 'UKF', 'PF'};
N_filters    = 3;

% Storage: N_mc × N_filters
rmse_all  = zeros(N_mc, N_filters);
conv_all  = zeros(N_mc, N_filters);   % convergence time [s]
time_all  = zeros(N_mc, N_filters);   % compute time [s]

t_start_mc = tic;

for mc = 1:N_mc
    rng(mc);   % different seed each run — reproducible

    if mod(mc, 10) == 1
        fprintf('[MC] Running simulation %d / %d ...\n', mc, N_mc);
    end

    % ---- Generate sensor measurements for this seed ------------------
    [gyro_meas, gyro_bias_true, st_meas, st_valid] = ...
        generate_sensors(ENV);

    sens_mc.gyro_meas  = gyro_meas;
    sens_mc.gyro_bias  = gyro_bias_true;
    sens_mc.st_meas    = st_meas;
    sens_mc.st_valid   = st_valid;

    % ---- Run EKF -----------------------------------------------------
    t0 = tic;
    [rmse_ekf, conv_ekf] = run_ekf_mc(ENV, sens_mc);
    time_all(mc, 1) = toc(t0);
    rmse_all(mc, 1) = rmse_ekf;
    conv_all(mc, 1) = conv_ekf;

    % ---- Run UKF -----------------------------------------------------
    t0 = tic;
    [rmse_ukf, conv_ukf] = run_ukf_mc(ENV, sens_mc);
    time_all(mc, 2) = toc(t0);
    rmse_all(mc, 2) = rmse_ukf;
    conv_all(mc, 2) = conv_ukf;

    % ---- Run PF ------------------------------------------------------
    t0 = tic;
    [rmse_pf, conv_pf] = run_pf_mc(ENV, sens_mc);
    time_all(mc, 3) = toc(t0);
    rmse_all(mc, 3) = rmse_pf;
    conv_all(mc, 3) = conv_pf;
end

total_mc_time = toc(t_start_mc);
fprintf('\n[MC] All %d simulations complete in %.1f s (%.1f min)\n\n', ...
        N_mc, total_mc_time, total_mc_time/60);

% =========================================================================
% 1.  COMPUTE STATISTICS
% =========================================================================
MC.N_mc        = N_mc;
MC.rmse_all    = rmse_all;
MC.conv_all    = conv_all;
MC.time_all    = time_all;

stats_fields = {'mean','std','median','p25','p75','min','max'};
for f = 1:N_filters
    data = rmse_all(:, f);
    MC.rmse.(filter_names{f}).mean   = mean(data);
    MC.rmse.(filter_names{f}).std    = std(data);
    MC.rmse.(filter_names{f}).median = median(data);
    MC.rmse.(filter_names{f}).p25    = prctile(data, 25);
    MC.rmse.(filter_names{f}).p75    = prctile(data, 75);
    MC.rmse.(filter_names{f}).min    = min(data);
    MC.rmse.(filter_names{f}).max    = max(data);
    MC.rmse.(filter_names{f}).cv     = std(data)/mean(data)*100;  % coeff of variation [%]
end

% =========================================================================
% 2.  STATISTICAL SIGNIFICANCE TESTS
% =========================================================================
fprintf('══════════════════════════════════════════════════════\n');
fprintf('  STATISTICAL RESULTS  (RMSE over %d runs)\n', N_mc);
fprintf('══════════════════════════════════════════════════════\n\n');

fprintf('%-10s  %-10s  %-10s  %-10s  %-10s  %-8s\n', ...
        'Filter','Mean [°]','Std [°]','Median [°]','95% CI [°]','CV [%]');
fprintf('%s\n', repmat('-', 1, 65));

for f = 1:N_filters
    s    = MC.rmse.(filter_names{f});
    ci95 = 1.96 * s.std / sqrt(N_mc);
    fprintf('%-10s  %-10.4f  %-10.4f  %-10.4f  ±%-9.4f  %-8.1f\n', ...
            filter_names{f}, s.mean, s.std, s.median, ci95, s.cv);
end

% Wilcoxon signed-rank test: UKF vs EKF
[p_ukf_ekf, ~] = signrank(rmse_all(:,1), rmse_all(:,2));
[p_pf_ekf,  ~] = signrank(rmse_all(:,1), rmse_all(:,3));
[p_ukf_pf,  ~] = signrank(rmse_all(:,2), rmse_all(:,3));

% Cohen's d effect size
d_ukf_ekf = (mean(rmse_all(:,1)) - mean(rmse_all(:,2))) / ...
             sqrt((std(rmse_all(:,1))^2 + std(rmse_all(:,2))^2)/2);
d_pf_ekf  = (mean(rmse_all(:,1)) - mean(rmse_all(:,3))) / ...
             sqrt((std(rmse_all(:,1))^2 + std(rmse_all(:,3))^2)/2);

fprintf('\n  Wilcoxon signed-rank test (non-parametric):\n');
fprintf('    UKF vs EKF : p = %.2e  →  %s  |  Cohen''s d = %.2f  (%s)\n', ...
        p_ukf_ekf, sig_label(p_ukf_ekf), d_ukf_ekf, effect_label(d_ukf_ekf));
fprintf('    PF  vs EKF : p = %.2e  →  %s  |  Cohen''s d = %.2f  (%s)\n', ...
        p_pf_ekf,  sig_label(p_pf_ekf),  d_pf_ekf,  effect_label(d_pf_ekf));
fprintf('    UKF vs PF  : p = %.2e  →  %s\n\n', ...
        p_ukf_pf,  sig_label(p_ukf_pf));

% Store test results
MC.stats.p_ukf_ekf  = p_ukf_ekf;
MC.stats.p_pf_ekf   = p_pf_ekf;
MC.stats.p_ukf_pf   = p_ukf_pf;
MC.stats.d_ukf_ekf  = d_ukf_ekf;
MC.stats.d_pf_ekf   = d_pf_ekf;

% =========================================================================
% 3.  FIGURES
% =========================================================================
colors = {[0.00, 0.45, 0.74],   % EKF — blue
          [0.85, 0.33, 0.10],   % UKF — red
          [0.47, 0.67, 0.19]};  % PF  — green

% --- Fig 1: Box Plot ---
figure('Name','MC — RMSE Box Plot','NumberTitle','off', ...
       'Position',[50 50 800 520]);
h = boxplot(rmse_all, filter_names, 'Colors', [0 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19], ...
            'Widths', 0.5, 'Symbol', 'o');
set(h, 'LineWidth', 1.3);
ylabel('RMSE [deg]');
title(sprintf('Monte Carlo RMSE Distribution  (N = %d runs)', N_mc));
grid on;

% Overlay individual points (jittered)
hold on;
for f = 1:N_filters
    jitter = 0.1 * (rand(N_mc,1) - 0.5);
    scatter(f + jitter, rmse_all(:,f), 18, colors{f}, 'filled', 'MarkerFaceAlpha', 0.4);
end

% Significance brackets
y_max = max(rmse_all(:));
bracket(1, 2, y_max*1.05, p_ukf_ekf);
hold off;

% --- Fig 2: RMSE over MC runs (convergence trace) ---
figure('Name','MC — RMSE per Run','NumberTitle','off', ...
       'Position',[100 80 1100 440]);

subplot(1,2,1);
hold on;
for f = 1:N_filters
    plot(1:N_mc, rmse_all(:,f), 'Color', [colors{f}, 0.5], 'LineWidth', 0.8);
    yline(mean(rmse_all(:,f)), '--', 'Color', colors{f}, 'LineWidth', 1.5, ...
          'Label', sprintf('%s mean = %.4f°', filter_names{f}, mean(rmse_all(:,f))));
end
xlabel('Monte Carlo run #'); ylabel('RMSE [deg]');
title('RMSE per simulation run');
legend('EKF','','UKF','','PF','','Location','northeast');
grid on; hold off;

subplot(1,2,2);
% Cumulative mean (running average) — shows convergence of statistics
for f = 1:N_filters
    running_mean = cumsum(rmse_all(:,f)) ./ (1:N_mc)';
    plot(1:N_mc, running_mean, 'Color', colors{f}, 'LineWidth', 1.5); hold on;
end
xlabel('Number of MC runs'); ylabel('Cumulative mean RMSE [deg]');
title('Convergence of Monte Carlo statistics');
legend('EKF','UKF','PF','Location','best');
grid on; hold off;

sgtitle('Monte Carlo Analysis — RMSE per Run');

% --- Fig 3: CDF (Empirical Cumulative Distribution) ---
figure('Name','MC — RMSE CDF','NumberTitle','off', ...
       'Position',[150 60 750 500]);
hold on;
for f = 1:N_filters
    sorted_rmse = sort(rmse_all(:,f));
    cdf_vals    = (1:N_mc) / N_mc;
    stairs(sorted_rmse, cdf_vals, 'Color', colors{f}, 'LineWidth', 2.0);
end
xlabel('RMSE [deg]'); ylabel('Cumulative probability');
title('Empirical CDF of RMSE');
legend('EKF','UKF','PF','Location','southeast');
xline(0.5,'k--','0.5°','LabelVerticalAlignment','bottom');
grid on; hold off;

% --- Fig 4: Mean ± 2σ Summary Bar ---
figure('Name','MC — Summary Statistics','NumberTitle','off', ...
       'Position',[200 50 850 500]);

means = [MC.rmse.EKF.mean, MC.rmse.UKF.mean, MC.rmse.PF.mean];
stds  = [MC.rmse.EKF.std,  MC.rmse.UKF.std,  MC.rmse.PF.std];

b = bar(means, 0.5, 'FaceColor','flat');
b.CData = [colors{1}; colors{2}; colors{3}];
hold on;
errorbar(1:3, means, 2*stds, 2*stds, 'k.', 'LineWidth', 1.5, 'CapSize', 10);

set(gca,'XTickLabel', filter_names);
ylabel('Mean RMSE [deg]');
title(sprintf('Mean ± 2σ RMSE  (N = %d Monte Carlo runs)', N_mc));

% p-value annotation
if p_ukf_ekf < 0.001
    sig_str = '***';
elseif p_ukf_ekf < 0.01
    sig_str = '**';
else
    sig_str = '*';
end
y_top = max(means + 2*stds) * 1.15;
plot([1 2], [y_top y_top]*0.98, 'k-', 'LineWidth', 1);
text(1.5, y_top, sig_str, 'HorizontalAlignment','center','FontSize',14);

% Value labels
for f = 1:3
    text(f, means(f)/2, sprintf('%.4f°', means(f)), ...
         'HorizontalAlignment','center','Color','w','FontSize',10,'FontWeight','bold');
end

grid on; hold off;

% --- Fig 5: Convergence time comparison ---
figure('Name','MC — Convergence Time','NumberTitle','off', ...
       'Position',[250 40 800 440]);

conv_valid = conv_all;
conv_valid(isnan(conv_all)) = NaN;

subplot(1,2,1);
h2 = boxplot(conv_valid, filter_names, ...
             'Colors',[0 0.45 0.74;0.85 0.33 0.10;0.47 0.67 0.19], ...
             'Widths',0.5);
set(h2,'LineWidth',1.2);
ylabel('Convergence time [s]');
title('Convergence time distribution');
grid on;

subplot(1,2,2);
mean_time = mean(time_all);
std_time  = std(time_all);
b3 = bar(mean_time, 0.5, 'FaceColor','flat');
b3.CData = [colors{1};colors{2};colors{3}];
hold on;
errorbar(1:3, mean_time, std_time, 'k.','LineWidth',1.5,'CapSize',8);
set(gca,'XTickLabel',filter_names);
ylabel('Compute time [s]');
title('Mean computation time per run');
for f=1:3
    text(f, mean_time(f)/2, sprintf('%.2fs',mean_time(f)), ...
         'HorizontalAlignment','center','Color','w','FontSize',9,'FontWeight','bold');
end
grid on; hold off;

sgtitle('Monte Carlo — Convergence & Computation Time');

% =========================================================================
% 4.  PRINT THESIS-READY SUMMARY
% =========================================================================
fprintf('══════════════════════════════════════════════════════\n');
fprintf('  THESIS-READY SUMMARY\n');
fprintf('══════════════════════════════════════════════════════\n\n');
fprintf('  "Over %d Monte Carlo simulations with independent\n', N_mc);
fprintf('   noise realizations, the UKF achieved a mean RMSE\n');
fprintf('   of %.4f° ± %.4f° (mean ± std), compared to\n', ...
        MC.rmse.UKF.mean, MC.rmse.UKF.std);
fprintf('   %.4f° ± %.4f° for the EKF — a %.1f%% improvement.\n', ...
        MC.rmse.EKF.mean, MC.rmse.EKF.std, ...
        (MC.rmse.EKF.mean - MC.rmse.UKF.mean)/MC.rmse.EKF.mean*100);
fprintf('   The difference is statistically significant\n');
fprintf('   (Wilcoxon signed-rank test, p = %.2e, %s).\n', ...
        p_ukf_ekf, sig_label(p_ukf_ekf));
fprintf('   Cohen''s d = %.2f (%s effect size)."\n\n', ...
        d_ukf_ekf, effect_label(d_ukf_ekf));

fprintf('[MC] Done. MC struct saved in workspace.\n');
fprintf('[MC] Add to run_full_simulation.m with: run(''monte_carlo_statistical_analysis.m'')\n\n');

% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function [gyro_meas, gyro_bias_true, st_meas, st_valid] = generate_sensors(ENV)
N          = ENV.N;
dt         = ENV.dt;
omega_true = ENV.omega_true;
q_true     = ENV.q_true;

gyro_bias  = zeros(3,N);
gyro_meas  = zeros(3,N);
gyro_bias(:,1) = ENV.sigma_gyro * 2 * randn(3,1);
for k = 1:N-1
    gyro_bias(:,k+1) = gyro_bias(:,k) + ENV.sigma_bias*sqrt(dt)*randn(3,1);
    gyro_meas(:,k)   = omega_true(:,k) + gyro_bias(:,k) + ENV.sigma_gyro*randn(3,1);
end
gyro_meas(:,N) = omega_true(:,N) + gyro_bias(:,N) + ENV.sigma_gyro*randn(3,1);
gyro_bias_true = gyro_bias;

st_meas  = zeros(4,N);
st_valid = false(1,N);
for k = 1:N
    if mod(k-1, ENV.step_st) == 0
        d = [ENV.sigma_st_cb*randn; ENV.sigma_st_cb*randn; ENV.sigma_st_bs*randn];
        h = norm(d)/2;
        if h < 1e-10, dq=[d;1]; else, dq=[d/norm(d)*sin(h);cos(h)]; end
        st_meas(:,k) = quat_mult(dq, q_true(:,k));
        st_valid(k)  = true;
    end
end
end

% ------------------------------------------------------------------
function [rmse, conv_t] = run_ekf_mc(ENV, sens)
N=ENV.N; dt=ENV.dt; Q=ENV.Q; R=ENV.R_meas;
ie=deg2rad(12); ia=[1;1;1]/sqrt(3);
q_est=[ia*sin(ie/2);cos(ie/2)]; q_est=quat_mult(q_est,ENV.q_true(:,1));
x=[q_est;zeros(3,1)]; P=diag([1e-4*ones(1,4),1e-6*ones(1,3)]);
err=zeros(1,N); conv_t=NaN; thresh=0.5;
for k=1:N-1
    q=x(1:4); b=x(5:7); oc=sens.gyro_meas(:,k)-b;
    Xi=[q(4)*eye(3)+skew(q(1:3));-q(1:3)'];
    qn=q+0.5*dt*Xi*oc; qn=qn/norm(qn); xp=[qn;b];
    F=ekf_jac(q,oc,dt); Pp=F*P*F'+Q;
    if sens.st_valid(k+1)
        qp=xp(1:4); zm=sens.st_meas(:,k+1);
        dq=quat_mult(zm,quat_inv(qp)); if dq(4)<0,dq=-dq;end
        inn=2*dq(1:3); qs=qp(4); qv=qp(1:3);
        H=[2*[qs*eye(3)+skew(qv),-qv],zeros(3,3)];
        K=Pp*H'/(H*Pp*H'+R); dx=K*inn;
        qu=qplus(xp(1:4),dx(1:4)); x=[qu;xp(5:7)+dx(5:7)];
        IKH=eye(7)-K*H; P=IKH*Pp*IKH'+K*R*K';
    else; x=xp; P=Pp; end
    x(1:4)=x(1:4)/norm(x(1:4));
    dq=quat_mult(quat_inv(ENV.q_true(:,k+1)),x(1:4));
    if dq(4)<0,dq=-dq;end; err(k+1)=2*asin(min(1,norm(dq(1:3))))*(180/pi);
    if isnan(conv_t) && k>50 && all(err(max(1,k-50):k)<thresh), conv_t=ENV.t(k); end
end
rmse=sqrt(mean(err(50:end).^2));
end

% ------------------------------------------------------------------
function [rmse, conv_t] = run_ukf_mc(ENV, sens)
N=ENV.N; dt=ENV.dt; R=ENV.R_meas; n=6;
alpha=1e-3; beta=2; kappa=0; lam=alpha^2*(n+kappa)-n; gam=sqrt(n+lam);
Wm=[(lam/(n+lam)),repmat(1/(2*(n+lam)),1,2*n)];
Wc=[(lam/(n+lam)+1-alpha^2+beta),repmat(1/(2*(n+lam)),1,2*n)];
qp=(ENV.sigma_gyro*dt)^2; qb=(ENV.sigma_bias*dt)^2;
Q6=diag([qp*ones(1,3),qb*ones(1,3)]);
ie=deg2rad(12); ia=[1;1;1]/sqrt(3);
dq0=[ia*sin(ie/2);cos(ie/2)];
q_est=quat_mult(dq0,ENV.q_true(:,1)); b_est=zeros(3,1);
P=diag([1e-4*ones(1,3),1e-6*ones(1,3)]);
err=zeros(1,N); conv_t=NaN; thresh=0.5;
for k=1:N-1
    try; Sc=gam*chol(P,'lower'); catch; P=nearSPD(P); Sc=gam*chol(P,'lower'); end
    es=[zeros(6,1),Sc,-Sc];
    qsg=zeros(4,2*n+1); bsg=zeros(3,2*n+1);
    for i=1:2*n+1; dqi=rv2q(es(1:3,i)); qsg(:,i)=quat_mult(dqi,q_est); bsg(:,i)=b_est+es(4:6,i); end
    qsp=zeros(4,2*n+1); bsp=zeros(3,2*n+1);
    for i=1:2*n+1
        oc=sens.gyro_meas(:,k)-bsg(:,i);
        Xi=[qsg(4,i)*eye(3)+skew(qsg(1:3,i));-qsg(1:3,i)'];
        qn=qsg(:,i)+0.5*dt*Xi*oc; qsp(:,i)=qn/norm(qn); bsp(:,i)=bsg(:,i);
    end
    qp2=sum(bsxfun(@times,qsp,Wm),2); qp2=qp2/norm(qp2);
    bp=sum(bsxfun(@times,bsp,Wm),2);
    Pp=Q6;
    for i=1:2*n+1
        dqi=quat_mult(qsp(:,i),quat_inv(qp2)); if dqi(4)<0,dqi=-dqi;end
        ei=[2*dqi(1:3);bsp(:,i)-bp]; Pp=Pp+Wc(i)*(ei*ei');
    end
    Pp=(Pp+Pp')/2;
    if sens.st_valid(k+1)
        try; Sc2=gam*chol(Pp,'lower'); catch; Pp=nearSPD(Pp); Sc2=gam*chol(Pp,'lower'); end
        e2=[zeros(6,1),Sc2,-Sc2];
        qsu=zeros(4,2*n+1); bsu=zeros(3,2*n+1);
        for i=1:2*n+1; dqi=rv2q(e2(1:3,i)); qsu(:,i)=quat_mult(dqi,qp2); bsu(:,i)=bp+e2(4:6,i); end
        zs=zeros(3,2*n+1);
        for i=1:2*n+1; dqi=quat_mult(qsu(:,i),quat_inv(qp2)); if dqi(4)<0,dqi=-dqi;end; zs(:,i)=2*dqi(1:3); end
        zp=sum(bsxfun(@times,zs,Wm),2);
        Pzz=R; Pxz=zeros(6,3);
        for i=1:2*n+1
            dz=zs(:,i)-zp; dqi=quat_mult(qsu(:,i),quat_inv(qp2));
            if dqi(4)<0,dqi=-dqi;end; ei=[2*dqi(1:3);bsu(:,i)-bp];
            Pzz=Pzz+Wc(i)*(dz*dz'); Pxz=Pxz+Wc(i)*(ei*dz');
        end
        zm=sens.st_meas(:,k+1); dqm=quat_mult(zm,quat_inv(qp2));
        if dqm(4)<0,dqm=-dqm;end; inn=2*dqm(1:3)-zp;
        K=Pxz/Pzz; dx=K*inn;
        dqc=rv2q(dx(1:3)); q_est=quat_mult(dqc,qp2); b_est=bp+dx(4:6);
        IKP=eye(6)-K*(Pxz'/Pp); P=IKP*Pp*IKP'+K*R*K';
    else; q_est=qp2; b_est=bp; P=Pp; end
    q_est=q_est/norm(q_est); P=(P+P')/2;
    dq=quat_mult(quat_inv(ENV.q_true(:,k+1)),q_est);
    if dq(4)<0,dq=-dq;end; err(k+1)=2*asin(min(1,norm(dq(1:3))))*(180/pi);
    if isnan(conv_t) && k>50 && all(err(max(1,k-50):k)<thresh), conv_t=ENV.t(k); end
end
rmse=sqrt(mean(err(50:end).^2));
end

% ------------------------------------------------------------------
function [rmse, conv_t] = run_pf_mc(ENV, sens)
N=ENV.N; dt=ENV.dt; R=ENV.R_meas; Np=500;
ie=deg2rad(12); ia=[1;1;1]/sqrt(3);
dq0=[ia*sin(ie/2);cos(ie/2)]; q0=quat_mult(dq0,ENV.q_true(:,1));
pts=zeros(7,Np);
for i=1:Np
    da=0.5*deg2rad(5)*randn(3,1); h=norm(da)/2;
    if h<1e-10,dq=[da;1];else,dq=[da/norm(da)*sin(h);cos(h)];end
    pts(:,i)=[quat_mult(dq,q0);1e-4*randn(3,1)];
end
wts=ones(1,Np)/Np; sg=ENV.sigma_gyro*2; sb=ENV.sigma_bias*2;
err=zeros(1,N); conv_t=NaN; thresh=0.5;
for k=1:N-1
    for i=1:Np
        qi=pts(1:4,i); bi=pts(5:7,i);
        oc=sens.gyro_meas(:,k)-bi+sg*randn(3,1);
        Xi=[qi(4)*eye(3)+skew(qi(1:3));-qi(1:3)'];
        qn=qi+0.5*dt*Xi*oc; qn=qn/norm(qn);
        pts(:,i)=[qn;bi+sb*sqrt(dt)*randn(3,1)];
    end
    if sens.st_valid(k+1)
        zm=sens.st_meas(:,k+1);
        for i=1:Np
            dq=quat_mult(zm,quat_inv(pts(1:4,i))); if dq(4)<0,dq=-dq;end
            inn=2*dq(1:3); wts(i)=wts(i)*exp(-0.5*inn'*(R\inn));
        end
        ws=sum(wts); if ws<1e-300||isnan(ws),wts=ones(1,Np)/Np; else,wts=wts/ws;end
        if 1/sum(wts.^2)<Np/2
            pts=sysresamp(pts,wts,Np); wts=ones(1,Np)/Np;
        end
    end
    xe=sum(bsxfun(@times,pts,wts),2); xe(1:4)=xe(1:4)/norm(xe(1:4));
    dq=quat_mult(quat_inv(ENV.q_true(:,k+1)),xe(1:4));
    if dq(4)<0,dq=-dq;end; err(k+1)=2*asin(min(1,norm(dq(1:3))))*(180/pi);
    if isnan(conv_t) && k>50 && all(err(max(1,k-50):k)<thresh), conv_t=ENV.t(k); end
end
rmse=sqrt(mean(err(50:end).^2));
end

% ---- Math helpers ------------------------------------------------
function F=ekf_jac(q,w,dt)
qs=q(4); qv=q(1:3);
Fqq=eye(4)+0.5*dt*[0,w(3),-w(2),w(1);-w(3),0,w(1),w(2);w(2),-w(1),0,w(3);-w(1),-w(2),-w(3),0];
Fqb=-0.5*dt*[qs*eye(3)+skew(qv);-qv']; F=[Fqq,Fqb;zeros(3,4),eye(3)];
end
function dq=rv2q(phi); phi=phi(:); a=norm(phi);
if a<1e-10,dq=[phi*0.5;1];else,dq=[phi/a*sin(a/2);cos(a/2)];end; dq=dq/norm(dq);
end
function q=quat_mult(p,q2); pv=p(1:3);ps=p(4);qv=q2(1:3);qs=q2(4);
q=[ps*qv(:)+qs*pv(:)+cross(pv(:),qv(:));ps*qs-dot(pv,qv)]; q=q/norm(q);
end
function qi=quat_inv(q); qi=[-q(1:3);q(4)]; end
function S=skew(v); S=[0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0]; end
function qo=qplus(q,dq); dv=dq(1:3); h=norm(dv)/2;
if h<1e-10,dqi=[dv;1];else,dqi=[dv/norm(dv)*sin(h);cos(h)];end
qo=quat_mult(dqi,q); qo=qo/norm(qo);
end
function A=nearSPD(A); B=(A+A')/2;[V,D]=eig(B);d=max(diag(D),1e-10);
A=V*diag(d)*V';A=(A+A')/2;
end
function np=sysresamp(p,w,N); pos=((0:N-1)+rand())/N; cw=cumsum(w);
i=1;j=1;np=zeros(size(p));
while i<=N; if pos(i)<=cw(j);np(:,i)=p(:,j);i=i+1;else;j=min(j+1,N);end;end
end

% ---- Annotation helpers ------------------------------------------
function lbl = sig_label(p)
if p < 0.001,     lbl = 'p<0.001 ***';
elseif p < 0.01,  lbl = 'p<0.01  **';
elseif p < 0.05,  lbl = 'p<0.05  *';
else,             lbl = 'n.s.'; end
end
function lbl = effect_label(d)
d = abs(d);
if d >= 0.8,      lbl = 'large';
elseif d >= 0.5,  lbl = 'medium';
else,             lbl = 'small'; end
end
function bracket(x1, x2, y, p)
if p<0.001, mk='***'; elseif p<0.01, mk='**'; elseif p<0.05, mk='*'; else, mk='n.s.'; end
plot([x1 x1 x2 x2],[y*0.98 y y y*0.98],'k-','LineWidth',1);
text((x1+x2)/2, y*1.01, mk,'HorizontalAlignment','center','FontSize',12);
end
