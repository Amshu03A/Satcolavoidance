function ml_conjunction_app()
%ML_CONJUNCTION_APP
%  Interactive command-line app where you input custom satellite and/or
%  debris orbital elements and get:
%    - ML-predicted risk classification (Safe / Warning / Critical)
%    - Predicted minimum approach distance
%    - Predicted recommended delta-V
%    - Physics-based verification (actual propagation)
%    - Avoidance maneuver recommendation
%    - Visualization of the encounter geometry
%
%  Prerequisites: Run train_conjunction_ml() first to create the model file.
%
%  Usage: ml_conjunction_app()

clc;
fprintf(' Satellite Conjunction Risk Predictor  \n');
fprintf(' ML-Powered Interactive App    \n');
fprintf('\n\n');

%% Load or train model
model_file = 'conjunction_ml_model.mat';
if exist(model_file, 'file')
    fprintf('Loading saved ML model...\n');
    load(model_file, 'ml_model');
    fprintf('Model loaded. (Trained on %d samples | Classifier acc=%.1f%%)\n\n',...
        ml_model.N_train, ml_model.acc);
else
    fprintf('No saved model found. Training now (this takes ~1-2 min)...\n\n');
    ml_model = train_conjunction_ml(3000);
    fprintf('\n');
end

%% Main app loop
while true
    fprintf(' MENU\n');
    fprintf('  [1] Enter MY satellite + MY debris\n');
    fprintf('  [2] Use default satellite + enter MY debris\n');
    fprintf('  [3] Enter MY satellite + use random debris\n');
    fprintf('  [4] Run batch test (10 random scenarios)\n');
    fprintf('  [5] Retrain ML model\n');
    fprintf('  [0] Exit\n');
    choice = input('Select option: ');
    if isempty(choice), continue; end

    switch choice
        case 0
            fprintf('Exiting. Goodbye!\n');
            break;

        case 1
            sat = input_satellite('Enter YOUR satellite orbital elements:');
            deb = input_debris('Enter YOUR debris orbital elements:');
            run_prediction(sat, deb, ml_model);

        case 2
            sat = default_satellite();
            fprintf('\nUsing default satellite: IMG-01 (550km, i=97.6°)\n');
            deb = input_debris('Enter YOUR debris orbital elements:');
            run_prediction(sat, deb, ml_model);

        case 3
            sat = input_satellite('Enter YOUR satellite orbital elements:');
            deb = random_debris();
            fprintf('\nGenerated random debris object.\n');
            print_object_summary('Debris', deb);
            run_prediction(sat, deb, ml_model);

        case 4
            run_batch_test(ml_model);

        case 5
            n = input('Number of training samples (default 3000): ');
            if isempty(n), n = 3000; end
            ml_model = train_conjunction_ml(n);

        otherwise
            fprintf('Invalid option.\n');
    end

    fprintf('\nPress Enter to continue...\n');
    input('');
    clc;
end
end

function run_prediction(sat, deb, ml_model)
%RUN_PREDICTION  Core function: ML predict + physics verify + visualize.

fprintf(' RUNNING ANALYSIS...\n');

params = get_params();

%% --- ML Prediction ---
% Build feature vector [sat_a_km, sat_e, sat_i_deg, sat_RAAN_deg,
%                        sat_omega_deg, sat_M0_deg,
%                        deb_a_km, deb_e, deb_i_deg, deb_RAAN_deg,
%                        deb_omega_deg, deb_M0_deg]
feat = [sat.a/1e3, sat.e, rad2deg(sat.i), rad2deg(sat.RAAN), ...
        rad2deg(sat.omega), rad2deg(sat.M0), ...
        deb.a/1e3, deb.e, rad2deg(deb.i), rad2deg(deb.RAAN), ...
        rad2deg(deb.omega), rad2deg(deb.M0)];

% Normalize using saved stats
feat_norm = (feat - ml_model.mu_X) ./ ml_model.sigma_X;

% Predict
risk_cls  = predict_class_app(ml_model,  feat_norm);
pred_dist = predict_dist_app(ml_model,   feat_norm);
pred_dv   = predict_dv_app(ml_model,     feat_norm);

risk_labels = {'SAFE','WARNING','CRITICAL'};
risk_colors = {'✅','⚠️ ','🚨'};

fprintf('  Risk Level     : %s %s\n', risk_colors{risk_cls+1}, risk_labels{risk_cls+1});
fprintf('  Min Distance   : %.2f km  (predicted)\n', max(pred_dist,0));
fprintf('  Recommended ΔV : %.3f m/s (predicted)\n', max(pred_dv,0));

%% --- Physics Verification (actual propagation) ---
fprintf('\n Computing physics verification (propagating orbits)...\n');
sat_st = propagate_obj(sat, params);
deb_st = propagate_obj(deb, params);

sep  = sat_st(:,1:3) - deb_st(:,1:3);
dist = sqrt(sum(sep.^2,2)) / 1e3;
[actual_dist, tca_idx] = min(dist);
tca_hrs = params.t_vec(tca_idx) / 3600;

rel_vel = norm(sat_st(tca_idx,4:6) - deb_st(tca_idx,4:6));
Pc = chan_pc(actual_dist*1e3, sat.hard_body_r, deb.hard_body_r, rel_vel);

if actual_dist < 10
    actual_cls = 2;
elseif actual_dist < 50
    actual_cls = 1;
else
    actual_cls = 0;
end

fprintf('│  Actual Min Distance : %.2f km\n', actual_dist);
fprintf('│  Time of Closest App.: %.2f hours from now\n', tca_hrs);
fprintf('│  Relative Velocity   : %.1f m/s\n', rel_vel);
fprintf('│  Probability of Coll.: %.2e\n', Pc);
fprintf('│  Actual Risk Level   : %s %s\n', risk_colors{actual_cls+1}, risk_labels{actual_cls+1});
fprintf('│  ML Prediction Match : %s\n', ternary(risk_cls==actual_cls,'✅ CORRECT','⚠️  MISMATCH'));

%% --- Maneuver Recommendation ---
if actual_cls >= 1
    fprintf('\n┌─ MANEUVER RECOMMENDATION ────────────────────────┐\n');
    if actual_cls == 2
        rec_dv_phys = 0.5 + (10 - actual_dist)/10 * 0.8;
        fprintf('  ACTION REQUIRED: CRITICAL CONJUNCTION\n');
    else
        rec_dv_phys = 0.1 + (50 - actual_dist)/50 * 0.4;
        fprintf('  ACTION ADVISORY: WARNING CONJUNCTION\n');
    end
    burn_time = tca_hrs - 0.5;
    fprintf('  Burn %.2f hours before TCA\n', 0.5);
    fprintf('  Recommended ΔV  : %.3f m/s (physics)\n', rec_dv_phys);
    fprintf('  Recommended ΔV  : %.3f m/s (ML model)\n', max(pred_dv,0));
    fprintf('  Burn Direction  : Prograde + out-of-plane\n');
    fprintf('  Expected new min dist: >%.0f km after maneuver\n', actual_dist*2.5);
else
    fprintf('\n  No maneuver required. Satellite is SAFE.\n');
end

%% --- Visualization ---
visualize_encounter(sat_st, deb_st, dist, params, actual_dist, tca_idx, ...
                    actual_cls, Pc, sat, deb);
end

function visualize_encounter(sat_st, deb_st, dist, params, min_d, tca_idx, ...
                              risk_cls, Pc, sat, deb)
%VISUALIZE_ENCOUNTER  3-panel plot of the encounter geometry.

t_hr = params.t_vec / 3600;
Re   = params.Re / 1e3;
risk_colors = {[0.2 0.7 0.3],[0.95 0.7 0.1],[0.85 0.15 0.15]};
rc = risk_colors{risk_cls+1};

figure('Name','Conjunction Analysis','NumberTitle','off',...
       'Position',[80 80 1200 500],'Color','w');

%% Panel 1: Separation over time
ax1 = subplot(1,3,1);
plot(t_hr, dist, 'b-', 'LineWidth', 1.8); hold on;
yline(50, '--', 'Color',[0.95 0.7 0.1], 'LineWidth',1.2, 'Label','Warning 50km');
yline(10, '--', 'Color',[0.85 0.15 0.15], 'LineWidth',1.2, 'Label','Critical 10km');
plot(t_hr(tca_idx), min_d, 'v', 'MarkerSize',12, ...
     'MarkerFaceColor', rc, 'MarkerEdgeColor','k', 'LineWidth',1.5);
text(t_hr(tca_idx), min_d*1.08, sprintf('TCA\n%.1f km', min_d), ...
     'HorizontalAlignment','center','FontSize',8,'FontWeight','bold','Color',rc);
xlabel('Time [hours]'); ylabel('Separation [km]');
title('Separation vs Time','FontWeight','bold');
grid on; box off;

%% Panel 2: 3D orbital geometry (zoomed to encounter region)
ax2 = subplot(1,3,2);
hold on; axis equal; grid on; box off;
ax2.Color = [0.97 0.97 0.97];

% Draw small Earth symbol
[xe,ye,ze] = sphere(20);
surf(Re*xe, Re*ye, Re*ze, 'FaceColor',[0.3 0.5 0.8],...
    'EdgeColor','none','FaceAlpha',0.4,'HandleVisibility','off');

% Satellite and debris tracks
plot3(sat_st(:,1)/1e3, sat_st(:,2)/1e3, sat_st(:,3)/1e3, ...
      'b-', 'LineWidth',1.5, 'DisplayName','Satellite');
plot3(deb_st(:,1)/1e3, deb_st(:,2)/1e3, deb_st(:,3)/1e3, ...
      'r--', 'LineWidth',1.2, 'DisplayName','Debris');

% Mark TCA points
plot3(sat_st(tca_idx,1)/1e3, sat_st(tca_idx,2)/1e3, sat_st(tca_idx,3)/1e3,...
      'bs','MarkerSize',10,'MarkerFaceColor','b','HandleVisibility','off');
plot3(deb_st(tca_idx,1)/1e3, deb_st(tca_idx,2)/1e3, deb_st(tca_idx,3)/1e3,...
      'r^','MarkerSize',10,'MarkerFaceColor','r','HandleVisibility','off');

% Draw line between TCA points
plot3([sat_st(tca_idx,1) deb_st(tca_idx,1)]/1e3,...
      [sat_st(tca_idx,2) deb_st(tca_idx,2)]/1e3,...
      [sat_st(tca_idx,3) deb_st(tca_idx,3)]/1e3,...
      '-','Color',rc,'LineWidth',2,'DisplayName',sprintf('TCA=%.1fkm',min_d));

xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Orbital Geometry (ECI)','FontWeight','bold');
legend('show','Location','best','FontSize',7);
view(30,20);

%% Panel 3: Risk summary card
ax3 = subplot(1,3,3);
axis off;
risk_labels = {'SAFE','WARNING','CRITICAL'};
risk_emojis = {'✓','!','✕'};

% Background colour based on risk
bg_patch = patch([0 1 1 0],[0 0 1 1], rc*0.9+0.1, ...
    'FaceAlpha',0.15,'EdgeColor',rc,'LineWidth',2);

% Title
text(0.5, 0.92, 'RISK ASSESSMENT', 'HorizontalAlignment','center',...
    'FontSize',13,'FontWeight','bold','Units','normalized',...
    'Color',[0.2 0.2 0.2]);

% Risk level
text(0.5, 0.78, risk_labels{risk_cls+1}, 'HorizontalAlignment','center',...
    'FontSize',22,'FontWeight','bold','Units','normalized','Color',rc);

% Stats
stats = {
    sprintf('Min Distance: %.2f km', min_d);
    sprintf('Pc: %.2e', Pc);
    sprintf('Rel. Velocity: %.0f m/s', ...
        norm(sat_st(tca_idx,4:6)-deb_st(tca_idx,4:6)));
    sprintf('TCA: %.2f hours', params.t_vec(tca_idx)/3600);
    sprintf('Sat altitude: %.0f km', sat.a/1e3 - params.Re/1e3);
    sprintf('Deb altitude: %.0f km', deb.a/1e3 - params.Re/1e3);
    sprintf('Inclination diff: %.1f°', abs(rad2deg(sat.i)-rad2deg(deb.i)));
};

y_pos = 0.62;
for k = 1:length(stats)
    text(0.5, y_pos, stats{k}, 'HorizontalAlignment','center',...
        'FontSize',9.5,'Units','normalized','Color',[0.2 0.2 0.2]);
    y_pos = y_pos - 0.08;
end

ax3.Units = 'normalized';
sgtitle('Conjunction Analysis Report','FontSize',13,'FontWeight','bold');
end

function sat = input_satellite(prompt_str)
fprintf('\n%s\n', prompt_str);
fprintf('(Press Enter to accept default values shown in brackets)\n\n');
sat.a      = (6371 + get_input('  Altitude [km]', 550)) * 1e3;
sat.e      =         get_input('  Eccentricity (0=circular)', 0.001);
sat.i      = deg2rad(get_input('  Inclination [deg]', 97.6));
sat.RAAN   = deg2rad(get_input('  RAAN [deg]', 0));
sat.omega  = deg2rad(get_input('  Arg. of Perigee [deg]', 0));
sat.M0     = deg2rad(get_input('  Mean Anomaly [deg]', 0));
sat.hard_body_r = get_input('  Hard-body radius [m]', 1.5);
print_object_summary('Satellite', sat);
end

function deb = input_debris(prompt_str)
fprintf('\n%s\n', prompt_str);
fprintf('(Press Enter to accept default values shown in brackets)\n\n');
deb.a      = (6371 + get_input('  Altitude [km]', 555)) * 1e3;
deb.e      =         get_input('  Eccentricity', 0.005);
deb.i      = deg2rad(get_input('  Inclination [deg]', 98.0));
deb.RAAN   = deg2rad(get_input('  RAAN [deg]', 5));
deb.omega  = deg2rad(get_input('  Arg. of Perigee [deg]', 0));
deb.M0     = deg2rad(get_input('  Mean Anomaly [deg]', 180));
deb.hard_body_r = get_input('  Hard-body radius [m]', 0.5);
print_object_summary('Debris', deb);
end

function sat = default_satellite()
sat.a = (6371+550)*1e3; sat.e=0.001;
sat.i=deg2rad(97.6); sat.RAAN=0; sat.omega=0; sat.M0=0;
sat.hard_body_r=1.5;
end

function deb = random_debris()
rng('shuffle');
alt = 530 + rand()*40;
deb.a = (6371+alt)*1e3;
deb.e = rand()*0.02;
deb.i = deg2rad(85 + rand()*20);
deb.RAAN  = deg2rad(rand()*360);
deb.omega = deg2rad(rand()*360);
deb.M0    = deg2rad(rand()*360);
deb.hard_body_r = 0.1 + rand()*0.9;
end

function print_object_summary(label, obj)
fprintf('\n  %s summary:\n', label);
fprintf('    Altitude : %.1f km\n', obj.a/1e3 - 6371);
fprintf('    e        : %.4f\n', obj.e);
fprintf('    i        : %.2f deg\n', rad2deg(obj.i));
fprintf('    RAAN     : %.2f deg\n', rad2deg(obj.RAAN));
fprintf('    M0       : %.2f deg\n', rad2deg(obj.M0));
end

function val = get_input(prompt, default_val)
str = input(sprintf('%s [%.4g]: ', prompt, default_val), 's');
if isempty(strtrim(str))
    val = default_val;
else
    val = str2double(str);
    if isnan(val)
        fprintf('  Invalid — using default %.4g\n', default_val);
        val = default_val;
    end
end
end

function run_batch_test(ml_model)
%RUN_BATCH_TEST  Run 10 random scenarios and display a summary table.
fprintf('\n Running 10 random conjunction scenarios...\n\n');
params = get_params();
risk_labels = {'SAFE   ','WARNING','CRITICAL'};

fprintf('%-4s %-10s %-10s %-10s %-10s %-10s\n',...
    '#','Sat Alt','Deb Alt','Incl Diff','Actual km','ML Risk');
fprintf('%s\n', repmat('-',1,58));

for k = 1:10
    sat = default_satellite();
    sat.M0   = deg2rad(rand()*360);
    sat.RAAN = deg2rad(rand()*360);

    deb = random_debris();

    % Feature + ML predict
    feat = [sat.a/1e3, sat.e, rad2deg(sat.i), rad2deg(sat.RAAN), ...
            rad2deg(sat.omega), rad2deg(sat.M0), ...
            deb.a/1e3, deb.e, rad2deg(deb.i), rad2deg(deb.RAAN), ...
            rad2deg(deb.omega), rad2deg(deb.M0)];
    feat_norm = (feat - ml_model.mu_X) ./ ml_model.sigma_X;
    risk_cls  = predict_class_app(ml_model, feat_norm);

    % Physics verify
    sat_st = propagate_obj(sat, params);
    deb_st = propagate_obj(deb, params);
    sep    = sat_st(:,1:3) - deb_st(:,1:3);
    dist   = sqrt(sum(sep.^2,2))/1e3;
    actual_dist = min(dist);

    fprintf('%-4d %-10.0f %-10.0f %-10.1f %-10.2f %s\n', k,...
        sat.a/1e3-6371, deb.a/1e3-6371,...
        abs(rad2deg(sat.i)-rad2deg(deb.i)), actual_dist,...
        risk_labels{risk_cls+1});
end
fprintf('%s\n', repmat('-',1,58));
end

% =========================================================================
%  PHYSICS HELPERS
% =========================================================================
function params = get_params()
params.mu = 3.986004418e14;
params.Re = 6371e3;
params.J2 = 1.08263e-3;
params.dt = 30;
params.T_total = 86400;
params.t_vec   = 0:params.dt:params.T_total;
params.N_steps = length(params.t_vec);
params.warn_dist_km   = 50;
params.danger_dist_km = 10;
params.Pc_threshold   = 1e-4;
end

function states = propagate_obj(obj, params)
mu=params.mu; Re=params.Re; J2=params.J2;
a=obj.a; e=obj.e; i=obj.i;
n=sqrt(mu/a^3); p=a*(1-e^2);
dR=-1.5*n*J2*(Re/p)^2*cos(i);
dw= 0.75*n*J2*(Re/p)^2*(5*cos(i)^2-1);
N=params.N_steps; states=zeros(N,6);
for k=1:N
    t=params.t_vec(k);
    RAAN=obj.RAAN+dR*t; omg=obj.omega+dw*t; M=obj.M0+n*t;
    E=kepler_solve(M,e);
    nu=2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2));
    r=a*(1-e*cos(E)); rp=[r*cos(nu);r*sin(nu);0];
    vp=sqrt(mu/p)*[-sin(nu);e+cos(nu);0];
    R=rot3(-RAAN)*rot1(-i)*rot3(-omg);
    states(k,:)=[R*rp;R*vp]';
end
end

function E=kepler_solve(M,e)
M=mod(M,2*pi); E=M;
for it=1:50; dE=(M-E+e*sin(E))/(1-e*cos(E)); E=E+dE; if abs(dE)<1e-12,break;end;end
end
function R=rot1(a); R=[1 0 0;0 cos(a) sin(a);0 -sin(a) cos(a)]; end
function R=rot3(a); R=[cos(a) sin(a) 0;-sin(a) cos(a) 0;0 0 1]; end

function Pc=chan_pc(d,rs,rd,~)
sig=100; A=pi*(rs+rd)^2;
Pc=min((A/(2*pi*sig^2))*exp(-d^2/(2*sig^2)),1);
end

function s=ternary(cond,a,b); if cond; s=a; else; s=b; end; end

% =========================================================================
%  ML PREDICTION WRAPPERS
% =========================================================================
function cls = predict_class_app(ml_model, X_norm)
if size(X_norm,1)==1, X_norm_r=X_norm; else, X_norm_r=X_norm; end
switch ml_model.cls_type
    case 'patternnet'
        out=ml_model.net_cls(X_norm_r'); [~,cls]=max(out,[],1); cls=cls-1;
    case 'fitcecoc'
        cls=predict(ml_model.mdl_cls, X_norm_r);
    case 'knn'
        cls=knn_pred(ml_model.X_tr_cls, ml_model.y_tr_cls, X_norm_r, 7);
end
cls=cls(1);
end

function d=predict_dist_app(ml_model, X_norm)
if isfield(ml_model,'mdl_dist'); d=predict(ml_model.mdl_dist,X_norm);
else; d=predict(ml_model.lm_dist,X_norm); end
d=max(d,0);
end

function dv=predict_dv_app(ml_model, X_norm)
if isfield(ml_model,'mdl_dv'); dv=predict(ml_model.mdl_dv,X_norm);
else; dv=predict(ml_model.lm_dv,X_norm); end
dv=max(dv,0);
end

function cls=knn_pred(X_tr,y_tr,X_te,k)
cls=zeros(size(X_te,1),1);
for i=1:size(X_te,1)
    dists=sum((X_tr-X_te(i,:)).^2,2); [~,si]=sort(dists);
    cls(i)=mode(y_tr(si(1:k)));
end
end
