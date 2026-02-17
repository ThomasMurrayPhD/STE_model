% Script to get model comparison metrics

clear;
close all

# % Use VBA toolbox for model comparison (ask Eddie)



%% load model fits
models(1).model_fits = importdata('model_PredictionBias\model_PredictionBias_fit.mat');
models(2).model_fits = importdata('model_PredictionBias\model_PredictionBias_fit2.mat');
% models(3).model_fits = importdata('model3_fit.mat');
% models(4).model_fits = importdata('model4_fit.mat');

model_names = {
    'Prediction bias (stim_noise)';
    'Prediction bias (PU)'};


N_models = numel(model_names);

%% Load excluded participants
fits = readtable('C:\Users\Tom\OneDrive - University of Cambridge\Cambridge\SpotTheEmotion_ThreatOfScream\DATA\fits.csv');
IDs = fits.Var1;
include = strcmp(fits.Include, 'True');
excluded_IDs = IDs(~include);

% remove IDs from model_fits
excluded_idx = ismember([models(1).model_fits.ID], excluded_IDs);
for i = 1:N_models
    models(i).model_fits(excluded_idx) = [];
end

N_inc = numel(models(1).model_fits);


%% Safe/Threat Concatenated
LMEs = nan(N_inc, N_models);
BICs = nan(N_inc, N_models);
AICs = nan(N_inc, N_models);
for iP = 1:N_inc
    for iM = 1:N_models
        if ~isempty(models(iM).model_fits(iP).est)
            LMEs(iP, iM) = models(iM).model_fits(iP).est.optim.LME;
            BICs(iP, iM) = models(iM).model_fits(iP).est.optim.BIC;
            AICs(iP, iM) = models(iM).model_fits(iP).est.optim.AIC;
        end
    end
end

valid = ~(isnan(LMEs) + isinf(LMEs));
LMEs_valid = LMEs(~any(~valid, 2), :);
BICs_valid = BICs(~any(~valid, 2), :);
AICs_valid = AICs(~any(~valid, 2), :);
BICs_sum = sum(BICs_valid, 1);
AICs_sum = sum(AICs_valid, 1);

[alpha,exp_r,xp,pxp,bor] = spm_BMS(LMEs_valid);

figure;
bar(1:N_models, pxp, 'EdgeColor', 'none');
set(gca, 'xticklabels', model_names, 'XTickLabelRotation', 45)
ylabel('PXP')
title('BMS')

fprintf('\nAll data:')
N_fit = N_inc - (sum(isnan(LMEs), 1) + sum(isinf(LMEs), 1));
fprintf('\n\tN fit:')
for i = 1:N_models
    fprintf('\n\t\t%s: %i/%i', model_names{i}, N_fit(i), N_inc);
end
fprintf('\n\tBayes Ombinus Risk: %1.3f', bor);
fprintf('\n\tAIC (sum):')
for i = 1:N_models
    fprintf('\n\t\t%s: %1.1f', model_names{i}, AICs_sum(i));
end
fprintf('\n\tBIC (sum):')
for i = 1:N_models
    fprintf('\n\t\t%s: %1.1f', model_names{i}, BICs_sum(i));
end


figure; hold on;
plot(1:N_models, LMEs_valid, 'color', [0 0.4470 0.7410, .2], 'linewidth', .1)
plot(1:N_models, mean(LMEs_valid, 1), '-o', 'color', 'black', 'linewidth', 2, 'MarkerFaceColor', [1,1,1])
set(gca, 'xtick', 1:N_models, 'xticklabels', model_names)
ylabel('LME')
title('All')

figure; hold on;
plot(1:N_models, BICs_valid, 'color', [0 0.4470 0.7410, .2], 'linewidth', .1)
plot(1:N_models, mean(BICs_valid, 1), '-o', 'color', 'black', 'linewidth', 2, 'MarkerFaceColor', [1,1,1])
set(gca, 'xtick', 1:N_models, 'xticklabels', model_names)
ylabel('BIC')
title('All')

figure; hold on;
plot(1:N_models, AICs_valid, 'color', [0 0.4470 0.7410, .2], 'linewidth', .1)
plot(1:N_models, mean(AICs_valid, 1), '-o', 'color', 'black', 'linewidth', 2, 'MarkerFaceColor', [1,1,1])
set(gca, 'xtick', 1:N_models, 'xticklabels', model_names)
ylabel('AIC')
title('All')

%% Safe only
safe_idx = find(arrayfun(@(x) strcmp(models(1).model_fits(x).condition, 'Safe'), 1:N_inc));
for i = 1:N_models
    safe_models(i).model_fits = models(i).model_fits(safe_idx);
end
N_safe = numel(safe_models(1).model_fits);

LMEs = nan(N_safe, N_models);
BICs = nan(N_safe, N_models);
AICs = nan(N_safe, N_models);
for iP = 1:N_safe
    for iM = 1:N_models
        if ~isempty(models(iM).model_fits(iP).est)
            LMEs(iP, iM) = safe_models(iM).model_fits(iP).est.optim.LME;
            BICs(iP, iM) = safe_models(iM).model_fits(iP).est.optim.BIC;
            AICs(iP, iM) = safe_models(iM).model_fits(iP).est.optim.AIC;
        end
    end
end

valid = ~(isnan(LMEs) + isinf(LMEs));
LMEs_valid = LMEs(~any(~valid, 2), :);
BICs_valid = BICs(~any(~valid, 2), :);
AICs_valid = AICs(~any(~valid, 2), :);
BICs_sum = sum(BICs_valid, 1);
AICs_sum = sum(AICs_valid, 1);
[alpha,exp_r,xp,pxp,bor] = spm_BMS(LMEs_valid);

figure;
bar(1:N_models, pxp, 'EdgeColor', 'none');
set(gca, 'xticklabels', model_names)
ylabel('PXP')
title('BMS (Safe only)')

fprintf('\n\nSafe only:')

N_fit = N_safe - (sum(isnan(LMEs), 1) + sum(isinf(LMEs), 1));
fprintf('\n\tN fit:')
for i = 1:N_models
    fprintf('\n\t\t%s: %i/%i', model_names{i}, N_fit(i), N_safe);
end
fprintf('\n\tBayes Ombinus Risk: %1.3f', bor);
fprintf('\n\tAIC (sum):')
for i = 1:N_models
    fprintf('\n\t\t%s: %1.1f', model_names{i}, AICs_sum(i));
end
fprintf('\n\tBIC (sum):')
for i = 1:N_models
    fprintf('\n\t\t%s: %1.1f', model_names{i}, BICs_sum(i));
end

figure; hold on;
plot(1:N_models, LMEs_valid, 'color', [0 0.4470 0.7410, .2], 'linewidth', .1)
plot(1:N_models, mean(LMEs_valid, 1), '-o', 'color', 'black', 'linewidth', 2, 'MarkerFaceColor', [1,1,1])
set(gca, 'xtick', 1:N_models, 'xticklabels', model_names)
ylabel('LME')
title('Safe only')


figure; hold on;
plot(1:N_models, BICs_valid, 'color', [0 0.4470 0.7410, .2], 'linewidth', .1)
plot(1:N_models, mean(BICs_valid, 1), '-o', 'color', 'black', 'linewidth', 2, 'MarkerFaceColor', [1,1,1])
set(gca, 'xtick', 1:N_models, 'xticklabels', model_names)
ylabel('BIC')
title('Safe only')

figure; hold on;
plot(1:N_models, AICs_valid, 'color', [0 0.4470 0.7410, .2], 'linewidth', .1)
plot(1:N_models, mean(AICs_valid, 1), '-o', 'color', 'black', 'linewidth', 2, 'MarkerFaceColor', [1,1,1])
set(gca, 'xtick', 1:N_models, 'xticklabels', model_names)
ylabel('AIC')
title('Safe only')


%% Threat only
threat_idx = find(arrayfun(@(x) strcmp(models(1).model_fits(x).condition, 'Threat'), 1:N_inc));
for i = 1:N_models
    threat_models(i).model_fits = models(i).model_fits(threat_idx);
end
N_threat = numel(threat_models(1).model_fits);

LMEs = nan(N_threat, N_models);
BICs = nan(N_threat, N_models);
AICs = nan(N_threat, N_models);
for iP = 1:N_threat
    for iM = 1:N_models
        if ~isempty(models(iM).model_fits(iP).est)
            LMEs(iP, iM) = threat_models(iM).model_fits(iP).est.optim.LME;
            BICs(iP, iM) = threat_models(iM).model_fits(iP).est.optim.BIC;
            AICs(iP, iM) = threat_models(iM).model_fits(iP).est.optim.AIC;
        end
    end
end

valid = ~(isnan(LMEs) + isinf(LMEs));
LMEs_valid = LMEs(~any(~valid, 2), :);
BICs_valid = BICs(~any(~valid, 2), :);
AICs_valid = AICs(~any(~valid, 2), :);
BICs_sum = sum(BICs_valid, 1);
AICs_sum = sum(AICs_valid, 1);
[alpha,exp_r,xp,pxp,bor] = spm_BMS(LMEs_valid);

figure;
bar(1:N_models, pxp, 'EdgeColor', 'none');
set(gca, 'xticklabels', model_names)
ylabel('PXP')
title('BMS (Threat only)')


fprintf('\n\nThreat only:')

N_fit = N_threat - (sum(isnan(LMEs), 1) + sum(isinf(LMEs), 1));
fprintf('\n\tN fit:')
for i = 1:N_models
    fprintf('\n\t\t%s: %i/%i', model_names{i}, N_fit(i), N_threat);
end
fprintf('\n\tBayes Ombinus Risk: %1.3f', bor);
fprintf('\n\tAIC (sum):')
for i = 1:N_models
    fprintf('\n\t\t%s: %1.1f', model_names{i}, AICs_sum(i));
end
fprintf('\n\tBIC (sum):')
for i = 1:N_models
    fprintf('\n\t\t%s: %1.1f', model_names{i}, BICs_sum(i));
end

figure; hold on;
plot(1:N_models, LMEs_valid, 'color', [0 0.4470 0.7410, .2], 'linewidth', .1)
plot(1:N_models, mean(LMEs_valid, 1), '-o', 'color', 'black', 'linewidth', 2, 'MarkerFaceColor', [1,1,1])
set(gca, 'xtick', 1:N_models, 'xticklabels', model_names)
ylabel('LME')
title('Threat only')

figure; hold on;
plot(1:N_models, BICs_valid, 'color', [0 0.4470 0.7410, .2], 'linewidth', .1)
plot(1:N_models, mean(BICs_valid, 1), '-o', 'color', 'black', 'linewidth', 2, 'MarkerFaceColor', [1,1,1])
set(gca, 'xtick', 1:N_models, 'xticklabels', model_names)
ylabel('BIC')
title('Threat only')

figure; hold on;
plot(1:N_models, AICs_valid, 'color', [0 0.4470 0.7410, .2], 'linewidth', .1)
plot(1:N_models, mean(AICs_valid, 1), '-o', 'color', 'black', 'linewidth', 2, 'MarkerFaceColor', [1,1,1])
set(gca, 'xtick', 1:N_models, 'xticklabels', model_names)
ylabel('AIC')
title('Threat only')

fprintf('\n\n')







%% Split obs models
model_fits = importdata('split_obs_fits.mat');
model_fits(excluded_idx) = []; % remove excluded IDs
N_inc = numel(model_fits);

% split model comparisons
LMEs_resp = nan(N_inc, 4);
LMEs_RT = nan(N_inc, 4);
for iP = 1:N_inc
    for iM = 1:4

        % response model
        if ~isempty(model_fits(iP).(['model', num2str(iM), '_resp_est']) )
            LMEs_resp(iP, iM) = model_fits(iP).(['model', num2str(iM), '_resp_est']).optim.LME;
        end
        
        % RT model
        if ~isempty(model_fits(iP).(['model', num2str(iM), '_RT_est']) )
            LMEs_RT(iP, iM) = model_fits(iP).(['model', num2str(iM), '_RT_est']).optim.LME;
        end

    end
end

% LMEs_resp = LMEs_resp(:, 1:3);
% LMEs_RT = LMEs_RT(:, 1:3);
% N_models = 3;
N_models = 4;

% response model BMS
valid_resp = ~(isnan(LMEs_resp) + isinf(LMEs_resp));
LMEs_resp_valid = LMEs_resp(~any(~valid_resp, 2), :);
[alpha,exp_r,xp,pxp,bor] = spm_BMS(LMEs_resp_valid);
figure('name', 'Response models (Split)');
bar(1:N_models, pxp, 'EdgeColor', 'none');
set(gca, 'xticklabels', model_names)
xtickangle(45)
ylabel('PXP')
title('BMS');

% RT model BMS
valid_RT = ~(isnan(LMEs_RT) + isinf(LMEs_RT));
LMEs_RT_valid = LMEs_RT(~any(~valid_RT, 2), :);
[alpha,exp_r,xp,pxp,bor] = spm_BMS(LMEs_RT_valid);
figure('name', 'RT models (Split)');
bar(1:N_models, pxp, 'EdgeColor', 'none');
set(gca, 'xticklabels', model_names)
xtickangle(45)
ylabel('PXP')
title('BMS');


