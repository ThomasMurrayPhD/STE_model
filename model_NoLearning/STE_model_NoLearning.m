% model 4 (no learning)


%%
clear;
close all;

addpath('..')


%% 

% example data (to get contingencies etc)
sub_data = readtable('..\STE_data\10369536_A_Threat.csv');

% input as %sad
u = sub_data.Outcome_p_sad/100;

% responses
y = [sub_data.Response_idx, log(sub_data.Response_RT)];


%% Get configuration structures
prc_model_config = tapas_ehgf_binary_config(); % perceptual model
obs_model_config = obs_NoLearning_comb_obs_config(); % response model
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;


%% Set parameters and priors

% fix prc params
prc_model_config.ommu(2) = -3;
prc_model_config.ommu(3) = -3;
prc_model_config.omsa(2) = 0;
prc_model_config.omsa(2) = 0;

% Obs model only
obs_model_config.logitalphamu = tapas_logit(.5, 1); % PSE
obs_model_config.logitalphasa = 4;

obs_model_config.logbetamu = log(15); % Slope
obs_model_config.logbetasa = 4;

obs_model_config.logitgamma0mu = tapas_logit(.0001, 1); % upper asymptote
obs_model_config.logitgamma0sa = 0;

obs_model_config.logitlambda0mu = tapas_logit(.0001, 1); % lower asymptote
obs_model_config.logitlambda0sa = 0;

obs_model_config.logitmumu = tapas_logit(.5, 1); % midpoint RT (Not used - uses alpha (PSE) instead)
obs_model_config.logitmusa = 0;

obs_model_config.logsigma0mu = log(.01); % width
obs_model_config.logsigma0sa = 4;

obs_model_config.loggamma1mu = log(0); % baseline logRT - fix to 0
obs_model_config.loggamma1sa = 0;

obs_model_config.loglambda1mu = log(7); % peak of logRT (above baseline)
obs_model_config.loglambda1sa = 4;

obs_model_config.logsigma1mu = log(1); % logRT noise
obs_model_config.logsigma1sa = 4;

obs_model_config = tapas_align_priors(obs_model_config);


%% simulate responses

r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = tapas_ehgf_binary_transp(r_temp, prc_model_config.priormus);

obs_params = obs_NoLearning_comb_obs_transp([], obs_model_config.priormus);

sim = tapas_simModel(u,...
    'tapas_ehgf_binary',...
    prc_params,...
    'obs_NoLearning_comb_obs',...
    obs_params,...
    123456789);

%% Visualise functions
% visualise_psychometric uses contingency space
N = 20;
all_y = nan(size(u,1), 2, N);
for i = 1:N
    sim = tapas_simModel(u,...
        'tapas_ehgf_binary',...
        prc_params,...
        'obs_NoLearning_comb_obs',...
        obs_params);
    all_y(:,:,i) = sim.y;
end

y_resp = squeeze(all_y(:,1,:)); % uxN
y_RT = squeeze(all_y(:,2,:));

figure;

% visualise psychometric
subplot(2,1,1); hold on;
all_resp_dists = zeros(N, 6);
for i = 1:N
    sim_resp_dist = arrayfun(@(x) mean(y_resp(sub_data.Outcome_p_sad==x, i)), 0:20:100);
    all_resp_dists(i, :) = sim_resp_dist;
    plot(0:20:100, sim_resp_dist, 'linewidth', 1, 'color', [.5,.5,.5]);
end
mean_resp = mean(all_resp_dists, 1);
plot(0:20:100, mean_resp, 'linewidth', 3, 'Color', [0    0.4470    0.7410]);
set(gca, 'Ylim', [0,1], 'Xtick', 0:20:100)
ylabel('p(Sad)')
title('simulated response distribution');

% visualise RT
subplot(2,1,2); hold on;
all_RT_dists = zeros(N, 6);
for i = 1:N
    sim_RT = y_RT(:, i);
    sim_RT_dist = arrayfun(@(x) mean(sim_RT(sub_data.Outcome_p_sad==x, 1)), 0:20:100);
    all_RT_dists(i, :) = sim_RT_dist;
    plot(0:20:100, sim_RT_dist, 'linewidth', 1, 'color', [.5,.5,.5]);
end
mean_RT = mean(all_RT_dists, 1);
plot(0:20:100, mean_RT, 'linewidth', 3, 'Color', [0    0.4470    0.7410]);


% G = @(x,mu,sigma0,gamma,lambda) gamma + lambda .* exp(-((x - mu).^2) ./ (2 * (sigma0.^2)));
% x=0:.01:1;
% gaus=G(x, tapas_sgm(obs_model_config.logitmumu, 1), exp(obs_model_config.logsigma0mu), exp(obs_model_config.loggamma1mu), exp(obs_model_config.loglambda1mu));
% plot(0:1:100, gaus, 'red', 'linewidth', 3);


set(gca, 'Xtick', 0:20:100)
xlabel('%Sad');
ylabel('logRT');
title('simulated RT distribution');

%% Fit
est = tapas_fitModel(...
    sim.y,...
    sim.u,...
    prc_model_config,...
    obs_model_config,...
    optim_config);


%% Parameter recovery
% So parameter recovery is pretty awful - keeps terminating optimisation.
% Let's try fixing gamma0 and lambda0 (upper/lower asymptotes of
% psychometric)
obs_model_config.logitgamma0mu = tapas_logit(.0001, 1);
obs_model_config.logitgamma0sa = 0;
obs_model_config.logitlambda0mu = tapas_logit(.0001, 1);
obs_model_config.logitlambda0sa = 0;
obs_model_config = tapas_align_priors(obs_model_config);

% set N iterations
N = 200;

% Parameters to recover
prc_param_names = {};
prc_param_idx   = [];
prc_param_space = {};

% obs_param_names = {'alpha', 'beta', 'gamma0', 'lambda0', 'mu', 'sigma0', 'gamma1', 'lambda1', 'sigma1'};
% obs_param_space = {'logit', 'log', 'logit', 'logit', 'logit', 'log', 'log', 'log', 'log'};
% obs_param_idx   = [1, 2, 3, 4, 5, 6, 7, 8, 9];
% obs_param_names = {'alpha', 'beta',  'mu', 'sigma0', 'gamma1', 'lambda1', 'sigma1'};
% obs_param_space = {'logit', 'log', 'logit', 'log', 'log', 'log', 'log'};
% obs_param_idx   = [1, 2, 5, 6, 7, 8, 9];

obs_param_names = {'alpha', 'beta', 'sigma0', 'lambda1', 'sigma1'};
obs_param_space = {'logit', 'log', 'log', 'log', 'log'};
obs_param_idx   = [1, 2, 6, 8, 9];


recov = parameter_recovery_master(u,...
    prc_model_config,...
    obs_model_config,...
    optim_config,...
    N,...
    prc_param_names,...
    prc_param_idx,...
    prc_param_space,...
    obs_param_names,...
    obs_param_idx,...
    obs_param_space);
save('model_NoLearning_recovery.mat', 'recov');

recovery_figures(recov);


%% Fit actual data

STE_dir = dir('..\STE_data\*.csv');
N_files = numel(STE_dir);
model_fits(N_files) = struct('ID', [], 'group', '', 'condition', '', 'est', struct());

% Just to be fancy
completion_times = zeros(N_files, 1);

for i = 1:N_files
    tic;
    fprintf('\nFitting dataset %i\n', i);
    if i>1
        avg_iter_time = mean(completion_times(1:i));
        fprintf('\tAverage iteration time = %1.2fs', avg_iter_time)
        estimated_total_time = avg_iter_time * ((N_files-i) + 1);
        fprintf('\n\tEstimated completion time = %im, %1.2fs\n\n', floor(estimated_total_time/60), rem(estimated_total_time,60));
    end

    % get ID, group, condition
    f_name = fullfile(STE_dir(i).folder, STE_dir(i).name);
    [~,name,~] = fileparts(f_name);
    tokens = split(name, '_');
    model_fits(i).ID           = str2double(tokens{1});
    model_fits(i).group        = tokens{2};
    model_fits(i).condition    = tokens{3};
    
    % load data
    sub_data = readtable(f_name);

    % model input
    u_sub = u;

    % subject responses
    sub_data.logRT = log(sub_data.Response_RT);
    y = [sub_data.Response_idx, sub_data.logRT]; % For learning-free model

    % remove missing
    missed = isnan(sub_data.Response_idx);
    u_sub = u_sub(~missed, :);
    y = y(~missed, :);

    % fit model
    try
        model_fits(i).est = tapas_fitModel(...
            y,...
            u_sub,...
            prc_model_config,...
            obs_model_config,...
            optim_config);
    catch
    end
    completion_times(i) = toc;
end

save('model_NoLearning_fit.mat', 'model_fits');





