function model = train_conjunction_ml(num_samples)
% Train three models on simulated satellite-debris data:
%   1. Classifier  - Safe / Warning / Critical risk
%   2. Regressor   - minimum approach distance [km]
%   3. Regressor   - recommended burn size [m/s]
%
% Tries Deep Learning Toolbox first (patternnet), falls back to
% Statistics Toolbox (fitcecoc / fitrensemble), then plain k-NN.
%
% Usage:
%   model = train_conjunction_ml(3000);
%   save('conjunction_ml_model.mat', 'ml_model');

if nargin < 1, num_samples = 3000; end

fprintf('==============================================\n');
fprintf('  Conjunction Risk ML Model Trainer\n');
fprintf('==============================================\n\n');

%% Generate data
[X, Y, ~] = generate_ml_training_data(num_samples);

min_dist = Y(:,1);
pc_vals  = Y(:,2);
risk_cls = Y(:,3);
rec_dv   = Y(:,4);

%% Normalize features
[X_norm, mu_feat, sig_feat] = zscore(X);

%% 80/20 train-test split
n      = size(X_norm, 1);
order  = randperm(n);
n_tr   = round(0.8*n);
tr     = order(1:n_tr);
te     = order(n_tr+1:end);

X_tr = X_norm(tr,:);  X_te = X_norm(te,:);
cls_tr = risk_cls(tr); cls_te = risk_cls(te);
dst_tr = min_dist(tr); dst_te = min_dist(te);
dv_tr  = rec_dv(tr);   dv_te  = rec_dv(te);

fprintf('\n[1/3] Training Risk Classifier (Safe / Warning / Critical)...\n');

has_deep  = (exist('patternnet','file') == 2);
has_stats = (exist('fitcecoc','file') == 2);

if has_deep
    fprintf('      Using Deep Learning Toolbox (patternnet)...\n');
    net = patternnet([64 32 16]);
    net.trainParam.showWindow    = false;
    net.trainParam.epochs        = 200;
    net.divideParam.trainRatio   = 0.85;
    net.divideParam.valRatio     = 0.15;
    net.divideParam.testRatio    = 0;
    t_hot = ind2vec(cls_tr'+1, 3);
    [net, ~] = train(net, X_tr', t_hot);
    model.cls_type = 'patternnet';
    model.net_cls  = net;
elseif has_stats
    fprintf('      Using Statistics Toolbox (fitcecoc + RBF SVM)...\n');
    tmpl = templateSVM('KernelFunction','rbf','Standardize',false);
    model.mdl_cls  = fitcecoc(X_tr, cls_tr, 'Learners', tmpl);
    model.cls_type = 'fitcecoc';
else
    fprintf('      Fallback: k-NN (no toolbox needed)...\n');
    model.X_tr_cls = X_tr;
    model.y_tr_cls = cls_tr;
    model.cls_type = 'knn';
end

preds_cls = run_classifier(model, X_te);
acc       = mean(preds_cls == cls_te) * 100;
fprintf('      Classifier accuracy: %.1f%%\n', acc);

% Confusion matrix with all three classes always present
C      = confusionmat(cls_te, preds_cls, 'Order', [0 1 2]);
labels = {'Safe','Warn','Crit'};
fprintf('      Confusion matrix (rows=true, cols=pred):\n');
fprintf('                Safe  Warn  Crit\n');
for r = 1:3
    fprintf('      %-6s  ', labels{r});
    fprintf('%5d ', C(r,:)); fprintf('\n');
end

%% Distance regressor
fprintf('\n[2/3] Training Min-Distance Regressor...\n');
if has_stats && exist('fitrensemble','file')
    model.mdl_dist = fitrensemble(X_tr, dst_tr, 'Method','LSBoost','NumLearningCycles',100);
    model.reg_type = 'ensemble';
else
    model.lm_dist  = fitlm(X_tr, dst_tr);
    model.reg_type = 'linear';
end
preds_dst = run_dist_reg(model, X_te);
rmse_dst  = sqrt(mean((preds_dst - dst_te).^2));
fprintf('      Distance RMSE: %.2f km\n', rmse_dst);

%% Delta-V regressor
fprintf('\n[3/3] Training Delta-V Regressor...\n');
if has_stats && exist('fitrensemble','file')
    model.mdl_dv = fitrensemble(X_tr, dv_tr, 'Method','LSBoost','NumLearningCycles',100);
else
    model.lm_dv  = fitlm(X_tr, dv_tr);
end
preds_dv = run_dv_reg(model, X_te);
rmse_dv  = sqrt(mean((preds_dv - dv_te).^2));
fprintf('      Delta-V RMSE: %.4f m/s\n', rmse_dv);

%% Store normalization and stats
model.mu_feat   = mu_feat;
model.sig_feat  = sig_feat;
model.acc       = acc;
model.rmse_dist = rmse_dst;
model.rmse_dv   = rmse_dv;
model.n_train   = n_tr;
model.n_test    = length(te);

ml_model = model;
save('conjunction_ml_model.mat', 'ml_model');
fprintf('\nModel saved to conjunction_ml_model.mat\n');

%% Diagnostic plots
figure('Name','ML Training Diagnostics','NumberTitle','off','Position',[100 100 1000 400]);

subplot(1,3,1);
confusionchart(cls_te, preds_cls, ...
    'RowSummary','row-normalized','ColumnSummary','column-normalized');
title(sprintf('Risk Classifier  (Acc=%.1f%%)', acc));

subplot(1,3,2);
scatter(dst_te, preds_dst, 20, 'filled','MarkerFaceAlpha',0.5);
hold on;
lims = [0, max([dst_te; preds_dst])*1.05];
plot(lims, lims, 'r--', 'LineWidth',1.5);
xlabel('True Min Distance [km]'); ylabel('Predicted [km]');
title(sprintf('Distance Regressor  RMSE=%.1f km', rmse_dst));
grid on;

subplot(1,3,3);
scatter(dv_te, preds_dv, 20, 'filled','MarkerFaceAlpha',0.5,...
    'MarkerFaceColor',[0.85 0.33 0.1]);
hold on;
lims2 = [0, max([dv_te; preds_dv])*1.05];
plot(lims2, lims2, 'r--', 'LineWidth',1.5);
xlabel('True \DeltaV [m/s]'); ylabel('Predicted [m/s]');
title(sprintf('\\DeltaV Regressor  RMSE=%.3f m/s', rmse_dv));
grid on;

sgtitle('ML Model Training Diagnostics','FontSize',13,'FontWeight','bold');

fprintf('\n==============================================\n');
fprintf('  Training Complete!\n');
fprintf('  Run ml_conjunction_app() to use the model.\n');
fprintf('==============================================\n');
end

% =========================================================================
%  Prediction helpers (also called by the app)
% =========================================================================
function cls = run_classifier(model, X)
switch model.cls_type
    case 'patternnet'
        out = model.net_cls(X');
        [~, cls] = max(out, [], 1);
        cls = cls' - 1;
    case 'fitcecoc'
        cls = predict(model.mdl_cls, X);
    case 'knn'
        cls = knn_vote(model.X_tr_cls, model.y_tr_cls, X, 7);
end
end

function d = run_dist_reg(model, X)
if isfield(model,'mdl_dist')
    d = predict(model.mdl_dist, X);
else
    d = predict(model.lm_dist, X);
end
end

function dv = run_dv_reg(model, X)
if isfield(model,'mdl_dv')
    dv = predict(model.mdl_dv, X);
else
    dv = predict(model.lm_dv, X);
end
dv = max(dv, 0);
end

function cls = knn_vote(X_tr, y_tr, X_te, k)
cls = zeros(size(X_te,1), 1);
for i = 1:size(X_te,1)
    dists      = sum((X_tr - X_te(i,:)).^2, 2);
    [~, order] = sort(dists);
    neighbors  = y_tr(order(1:k));
    cls(i)     = mode(neighbors);
end
end
