% =========================================================================
% stress_tests_and_robustness.m
% ADCS Filter Comparison Project — Module 5: Stress Tests
%
% Runs all 5 scenarios (S1–S5) and compares filter robustness.
% Each scenario re-runs M1 sensor sim with modified noise, then
% re-runs all 3 filters and records RMSE.
%
% Scenarios:
%   S1 — Nominal     (baseline)
%   S2 — High IMU    (gyro noise ×10)
%   S3 — Degraded ST (star tracker noise ×6)
%   S4 — Both high   (IMU ×10 + ST ×6)
%   S5 — ST outage   (blind 10s every 60s)
%
% Output: STRESS struct + summary heatmap
%
% Author  : Hatem Adel Saber
% =========================================================================

if ~exist('ENV','var')
    error('[M5] Run M0 first to generate ENV.');
end

fprintf('\n========== MODULE 5: STRESS TESTS ==========\n\n');

scenarios = {'S1 Nominal','S2 High IMU','S3 Degraded ST','S4 Both high','S5 ST outage'};
N_scen    = length(scenarios);
filter_names = {'EKF','UKF','PF'};
rmse_table   = zeros(N_scen, 3);   % rows=scenarios, cols=filters

for s = 1:N_scen
    fprintf('--- Running %s ---\n', scenarios{s});

    % Build scenario-specific ENV copy
    env_s = ENV;

    switch s
        case 1   % S1 Nominal — no change
            % unchanged

        case 2   % S2 High IMU noise ×10
            env_s.sigma_gyro = ENV.sigma_gyro * ENV.stress.high_imu_mult;
            env_s.sigma_bias = ENV.sigma_bias * ENV.stress.high_imu_mult;
            q_q = (env_s.sigma_gyro * env_s.dt)^2;
            q_b = (env_s.sigma_bias * env_s.dt)^2;
            env_s.Q = diag([q_q*ones(1,4), q_b*ones(1,3)]);

        case 3   % S3 Degraded Star Tracker ×6
            env_s.sigma_st_cb = ENV.sigma_st_cb * ENV.stress.high_st_mult;
            env_s.sigma_st_bs = ENV.sigma_st_bs * ENV.stress.high_st_mult;
            env_s.R_meas = diag([env_s.sigma_st_cb^2, ...
                                  env_s.sigma_st_cb^2, ...
                                  env_s.sigma_st_bs^2]);

        case 4   % S4 Both high
            env_s.sigma_gyro  = ENV.sigma_gyro * ENV.stress.high_imu_mult;
            env_s.sigma_bias  = ENV.sigma_bias * ENV.stress.high_imu_mult;
            env_s.sigma_st_cb = ENV.sigma_st_cb * ENV.stress.high_st_mult;
            env_s.sigma_st_bs = ENV.sigma_st_bs * ENV.stress.high_st_mult;
            q_q = (env_s.sigma_gyro*env_s.dt)^2;
            q_b = (env_s.sigma_bias*env_s.dt)^2;
            env_s.Q = diag([q_q*ones(1,4), q_b*ones(1,3)]);
            env_s.R_meas = diag([env_s.sigma_st_cb^2, ...
                                  env_s.sigma_st_cb^2, ...
                                  env_s.sigma_st_bs^2]);

        case 5   % S5 ST outage — nominal noise but periodic blind spots
            % handled in sensor sim below
    end

    % ---- Generate sensor measurements for this scenario ---------------
    N   = env_s.N;
    dt  = env_s.dt;
    q_true     = env_s.q_true;
    omega_true = env_s.omega_true;

    % IMU
    gyro_bias = zeros(3,N);
    gyro_meas = zeros(3,N);
    gyro_bias(:,1) = env_s.sigma_gyro * 2 * randn(3,1);
    for k = 1:N-1
        gyro_bias(:,k+1) = gyro_bias(:,k) + env_s.sigma_bias*sqrt(dt)*randn(3,1);
        gyro_meas(:,k)   = omega_true(:,k) + gyro_bias(:,k) + env_s.sigma_gyro*randn(3,1);
    end
    gyro_meas(:,N) = omega_true(:,N) + gyro_bias(:,N) + env_s.sigma_gyro*randn(3,1);

    % Star Tracker
    st_meas  = zeros(4,N);
    st_valid = false(1,N);
    for k = 1:N
        is_outage = false;
        if s == 5
            % blind for outage_dur steps every outage_period steps
            phase = mod(k-1, env_s.stress.outage_period);
            if phase < env_s.stress.outage_dur
                is_outage = true;
            end
        end

        if ~is_outage && mod(k-1, env_s.step_st) == 0
            d_angle = [env_s.sigma_st_cb*randn; ...
                       env_s.sigma_st_cb*randn; ...
                       env_s.sigma_st_bs*randn];
            half = norm(d_angle)/2;
            if half < 1e-10, dq = [d_angle;1];
            else, dq = [d_angle/norm(d_angle)*sin(half); cos(half)]; end
            st_meas(:,k) = quat_mult(dq, q_true(:,k));
            st_valid(k)  = true;
        end
    end

    sens_s.gyro_meas = gyro_meas;
    sens_s.gyro_bias = gyro_bias;
    sens_s.st_meas   = st_meas;
    sens_s.st_valid  = st_valid;

    % ---- Run EKF -------------------------------------------------------
    rmse_table(s,1) = run_ekf_stress(env_s, sens_s);
    fprintf('  EKF RMSE = %.4f deg\n', rmse_table(s,1));

    % ---- Run UKF -------------------------------------------------------
    rmse_table(s,2) = run_ukf_stress(env_s, sens_s);
    fprintf('  UKF RMSE = %.4f deg\n', rmse_table(s,2));

    % ---- Run PF --------------------------------------------------------
    rmse_table(s,3) = run_pf_stress(env_s, sens_s);
    fprintf('  PF  RMSE = %.4f deg\n\n', rmse_table(s,3));
end

% -------------------------------------------------------------------------
% RESULTS TABLE
% -------------------------------------------------------------------------
fprintf('\n========== STRESS TEST RESULTS (RMSE deg) ==========\n');
fprintf('%-20s  %-10s  %-10s  %-10s\n','Scenario','EKF','UKF','PF');
fprintf('%s\n', repmat('-',1,54));
for s = 1:N_scen
    fprintf('%-20s  %-10.4f  %-10.4f  %-10.4f\n', ...
            scenarios{s}, rmse_table(s,1), rmse_table(s,2), rmse_table(s,3));
end

STRESS.scenarios  = scenarios;
STRESS.rmse_table = rmse_table;

% -------------------------------------------------------------------------
% HEATMAP VISUALIZATION
% -------------------------------------------------------------------------
figure('Name','M5 — Stress Test Heatmap','NumberTitle','off', ...
       'Position',[100 100 700 450]);

imagesc(rmse_table);
colormap(flipud(hot));
colorbar;
set(gca, 'XTickLabel', filter_names, 'XTick', 1:3, ...
         'YTickLabel', scenarios,    'YTick', 1:N_scen);
xlabel('Filter'); ylabel('Scenario');
title('RMSE [deg] — Filter × Scenario (darker = worse)');

for s = 1:N_scen
    for f = 1:3
        text(f, s, sprintf('%.4f', rmse_table(s,f)), ...
             'HorizontalAlignment','center','FontSize',9,'Color','w','FontWeight','bold');
    end
end

% -------------------------------------------------------------------------
% BAR GROUP COMPARISON
% -------------------------------------------------------------------------
figure('Name','M5 — RMSE Stress Comparison','NumberTitle','off', ...
       'Position',[120 80 1000 420]);

bar(rmse_table, 0.7);
set(gca, 'XTickLabel', scenarios);
xtickangle(15);
ylabel('RMSE [deg]');
title('Filter robustness across all stress scenarios');
legend('EKF','UKF','PF','Location','northwest');
grid on;

fprintf('\n[M5] Stress tests complete. STRESS struct ready.\n');
fprintf('[M5] Key insight: check S4 and S5 — these reveal where EKF breaks.\n\n');

% =========================================================================
% INLINE FILTER RUNNERS  (self-contained, use env_s and sens_s)
% =========================================================================

function rmse = run_ekf_stress(env_s, sens_s)
N=env_s.N; dt=env_s.dt; Q=env_s.Q; R_meas=env_s.R_meas;
ie=deg2rad(12); ia=[1;1;1]/norm([1;1;1]);
dq0=[ia*sin(ie/2);cos(ie/2)];
q0=quat_mult(dq0,env_s.q_true(:,1));
x_est=[q0;zeros(3,1)];
P=diag([1e-4*ones(1,4),1e-6*ones(1,3)]);
err=zeros(1,N);
for k=1:N-1
    q=x_est(1:4); b=x_est(5:7);
    oc=sens_s.gyro_meas(:,k)-b;
    Xi=[q(4)*eye(3)+skew(q(1:3));-q(1:3)'];
    qn=q+0.5*dt*Xi*oc; qn=qn/norm(qn);
    xp=[qn;b];
    F=ekf_F(q,oc,dt);
    Pp=F*P*F'+Q;
    if sens_s.st_valid(k+1)
        qp=xp(1:4); zm=sens_s.st_meas(:,k+1);
        dq=quat_mult(zm,quat_inv(qp)); if dq(4)<0,dq=-dq;end
        inn=2*dq(1:3);
        qs=qp(4); qv=qp(1:3);
        Hq=2*[qs*eye(3)+skew(qv),-qv]; H=[Hq,zeros(3,3)];
        Si=H*Pp*H'+R_meas; K=Pp*H'/Si;
        dx=K*inn; qupd=qplus(xp(1:4),dx(1:4)); bupd=xp(5:7)+dx(5:7);
        x_est=[qupd;bupd]; IKH=eye(7)-K*H; P=IKH*Pp*IKH'+K*R_meas*K';
    else; x_est=xp; P=Pp; end
    x_est(1:4)=x_est(1:4)/norm(x_est(1:4));
    dq=quat_mult(quat_inv(env_s.q_true(:,k+1)),x_est(1:4));
    if dq(4)<0,dq=-dq;end
    err(k+1)=2*asin(min(1,norm(dq(1:3))))*(180/pi);
end
rmse=sqrt(mean(err(50:end).^2));
end

function rmse = run_ukf_stress(env_s, sens_s)
% UKF stress runner — 6-state error formulation (matches unscented_kalman_filter.m)
N=env_s.N; dt=env_s.dt; R_meas=env_s.R_meas;
n=6;  % error-state: 3 rotation + 3 bias
alpha=1e-3; beta=2; kappa=0;
lam=alpha^2*(n+kappa)-n; gam=sqrt(n+lam);
Wm=[(lam/(n+lam)), repmat(1/(2*(n+lam)),1,2*n)];
Wc=[(lam/(n+lam)+1-alpha^2+beta), repmat(1/(2*(n+lam)),1,2*n)];

% Process noise 6×6
q_phi=(env_s.sigma_gyro*dt)^2; q_b=(env_s.sigma_bias*dt)^2;
Q6=diag([q_phi*ones(1,3), q_b*ones(1,3)]);

ie=deg2rad(12); ia=[1;1;1]/norm([1;1;1]);
dq0=[ia*sin(ie/2);cos(ie/2)];
q_est=quat_mult(dq0,env_s.q_true(:,1)); b_est=zeros(3,1);
P=diag([1e-4*ones(1,3),1e-6*ones(1,3)]);  % 6×6
err=zeros(1,N);

for k=1:N-1
    % --- PREDICT ---
    try; Sc=gam*chol(P,'lower'); catch; P=nearSPD(P); Sc=gam*chol(P,'lower'); end
    eps_s=[zeros(6,1), Sc, -Sc];  % 6×(2n+1)
    q_sg=zeros(4,2*n+1); b_sg=zeros(3,2*n+1);
    for i=1:2*n+1
        dqi=rv2q(eps_s(1:3,i)); q_sg(:,i)=quat_mult(dqi,q_est); b_sg(:,i)=b_est+eps_s(4:6,i);
    end
    qsp=zeros(4,2*n+1); bsp=zeros(3,2*n+1);
    for i=1:2*n+1
        oc=sens_s.gyro_meas(:,k)-b_sg(:,i);
        Xi=[q_sg(4,i)*eye(3)+skew(q_sg(1:3,i));-q_sg(1:3,i)'];
        qn=q_sg(:,i)+0.5*dt*Xi*oc; qsp(:,i)=qn/norm(qn); bsp(:,i)=b_sg(:,i);
    end
    qp=sum(bsxfun(@times,qsp,Wm),2); qp=qp/norm(qp);
    bp=sum(bsxfun(@times,bsp,Wm),2);
    Pp=Q6;
    for i=1:2*n+1
        dqi=quat_mult(qsp(:,i),quat_inv(qp)); if dqi(4)<0,dqi=-dqi;end
        ei=[2*dqi(1:3);bsp(:,i)-bp]; Pp=Pp+Wc(i)*(ei*ei');
    end
    Pp=(Pp+Pp')/2;

    % --- UPDATE ---
    if sens_s.st_valid(k+1)
        try; Sc2=gam*chol(Pp,'lower'); catch; Pp=nearSPD(Pp); Sc2=gam*chol(Pp,'lower'); end
        eps2=[zeros(6,1),Sc2,-Sc2];
        qsu=zeros(4,2*n+1); bsu=zeros(3,2*n+1);
        for i=1:2*n+1
            dqi=rv2q(eps2(1:3,i)); qsu(:,i)=quat_mult(dqi,qp); bsu(:,i)=bp+eps2(4:6,i);
        end
        zs=zeros(3,2*n+1);
        for i=1:2*n+1
            dqi=quat_mult(qsu(:,i),quat_inv(qp)); if dqi(4)<0,dqi=-dqi;end; zs(:,i)=2*dqi(1:3);
        end
        zp=sum(bsxfun(@times,zs,Wm),2);
        Pzz=R_meas; Pxz=zeros(6,3);
        for i=1:2*n+1
            dz=zs(:,i)-zp; dqi=quat_mult(qsu(:,i),quat_inv(qp));
            if dqi(4)<0,dqi=-dqi;end
            ei=[2*dqi(1:3);bsu(:,i)-bp];
            Pzz=Pzz+Wc(i)*(dz*dz'); Pxz=Pxz+Wc(i)*(ei*dz');
        end
        zm=sens_s.st_meas(:,k+1);
        dqm=quat_mult(zm,quat_inv(qp)); if dqm(4)<0,dqm=-dqm;end
        inn=2*dqm(1:3)-zp;
        K=Pxz/Pzz; dx=K*inn;  % 6×1
        dqc=rv2q(dx(1:3)); q_est=quat_mult(dqc,qp); b_est=bp+dx(4:6);
        IKP=eye(6)-K*(Pxz'/Pp); P=IKP*Pp*IKP'+K*R_meas*K';
    else
        q_est=qp; b_est=bp; P=Pp;
    end
    q_est=q_est/norm(q_est); P=(P+P')/2;
    dq=quat_mult(quat_inv(env_s.q_true(:,k+1)),q_est);
    if dq(4)<0,dq=-dq;end; err(k+1)=2*asin(min(1,norm(dq(1:3))))*(180/pi);
end
rmse=sqrt(mean(err(50:end).^2));
end

function rmse = run_pf_stress(env_s, sens_s)
N=env_s.N; dt=env_s.dt; R_meas=env_s.R_meas; Np=300;
ie=deg2rad(12); ia=[1;1;1]/norm([1;1;1]);
dq0=[ia*sin(ie/2);cos(ie/2)]; q0=quat_mult(dq0,env_s.q_true(:,1));
parts=zeros(7,Np);
for i=1:Np
    da=0.5*deg2rad(5)*randn(3,1); hf=norm(da)/2;
    if hf<1e-10,dq=[da;1];else,dq=[da/norm(da)*sin(hf);cos(hf)];end
    parts(:,i)=[quat_mult(dq,q0);1e-4*randn(3,1)];
end
wts=ones(1,Np)/Np; sg=env_s.sigma_gyro*2; sb=env_s.sigma_bias*2;
err=zeros(1,N);
for k=1:N-1
    for i=1:Np
        qi=parts(1:4,i); bi=parts(5:7,i);
        oc=sens_s.gyro_meas(:,k)-bi+sg*randn(3,1);
        Xi=[qi(4)*eye(3)+skew(qi(1:3));-qi(1:3)'];
        qn=qi+0.5*dt*Xi*oc; qn=qn/norm(qn);
        parts(:,i)=[qn;bi+sb*sqrt(dt)*randn(3,1)];
    end
    if sens_s.st_valid(k+1)
        zm=sens_s.st_meas(:,k+1);
        xm=sum(bsxfun(@times,parts,wts),2); xm(1:4)=xm(1:4)/norm(xm(1:4));
        for i=1:Np
            dq=quat_mult(zm,quat_inv(parts(1:4,i))); if dq(4)<0,dq=-dq;end
            inn=2*dq(1:3); lw=-0.5*inn'*(R_meas\inn); wts(i)=wts(i)*exp(lw);
        end
        ws=sum(wts); if ws<1e-300||isnan(ws),wts=ones(1,Np)/Np; else,wts=wts/ws;end
        if 1/sum(wts.^2)<Np/2
            parts=sysresamp(parts,wts,Np); wts=ones(1,Np)/Np;
        end
    end
    xe=sum(bsxfun(@times,parts,wts),2); xe(1:4)=xe(1:4)/norm(xe(1:4));
    dq=quat_mult(quat_inv(env_s.q_true(:,k+1)),xe(1:4));
    if dq(4)<0,dq=-dq;end; err(k+1)=2*asin(min(1,norm(dq(1:3))))*(180/pi);
end
rmse=sqrt(mean(err(50:end).^2));
end

% ---- Shared helpers ----
function F=ekf_F(q,w,dt)
qs=q(4); qv=q(1:3);
Fqq=eye(4)+0.5*dt*[0,w(3),-w(2),w(1);-w(3),0,w(1),w(2);w(2),-w(1),0,w(3);-w(1),-w(2),-w(3),0];
Fqb=-0.5*dt*[qs*eye(3)+skew(qv);-qv']; F=[Fqq,Fqb;zeros(3,4),eye(3)];
end
function dq=rv2q(phi); phi=phi(:); a=norm(phi);
if a<1e-10, dq=[phi*0.5;1]; else, dq=[phi/a*sin(a/2);cos(a/2)]; end
dq=dq/norm(dq); end
function q=quat_mult(p,q2); pv=p(1:3);ps=p(4);qv=q2(1:3);qs=q2(4);
q=[ps*qv(:)+qs*pv(:)+cross(pv(:),qv(:));ps*qs-dot(pv,qv)]; q=q/norm(q); end
function qi=quat_inv(q); qi=[-q(1:3);q(4)]; end
function S=skew(v); S=[0,-v(3),v(2);v(3),0,-v(1);-v(2),v(1),0]; end
function qo=qplus(q,dq); dv=dq(1:3); h=norm(dv)/2;
if h<1e-10,dqi=[dv;1];else,dqi=[dv/norm(dv)*sin(h);cos(h)];end
qo=quat_mult(dqi,q); qo=qo/norm(qo); end
function A=nearSPD(A); B=(A+A')/2; [V,D]=eig(B); d=max(diag(D),1e-10);
A=V*diag(d)*V'; A=(A+A')/2; end
function np=sysresamp(p,w,N); pos=((0:N-1)+rand())/N; cw=cumsum(w);
i=1;j=1;np=zeros(size(p));
while i<=N; if pos(i)<=cw(j);np(:,i)=p(:,j);i=i+1;else;j=min(j+1,N);end;end;end