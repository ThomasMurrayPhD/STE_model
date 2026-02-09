% script to play around with recovery of noise parameter in logRT model

close all;

% example data (to get contingencies etc)
sub_data = readtable('..\STE_data\10369536_A_Threat.csv');

% input
u = double(sub_data.Cue_idx == sub_data.Outcome_idx);

% model configs
prc_model_config = tapas_hgf_binary_config;
obs_model_config = tapas_logrt_linear_binary_config;
optim_config = tapas_quasinewton_optim_config;
optim_config.nRandInit = 1;

% 2 level
prc_model_config.logkamu(2) = log(0);
prc_model_config.logkasa(2) = 0;
prc_model_config.omsa(3) = 0;
prc_model_config = tapas_align_priors(prc_model_config);

% get parameters
r = [];
r.c_prc.n_levels = 3;
prc_params = tapas_hgf_binary_transp(r, prc_model_config.priormus);
obs_params = tapas_logrt_linear_binary_transp([], obs_model_config.priormus);

% simulate
sim = tapas_simModel(u, 'tapas_hgf_binary', prc_params, 'tapas_logrt_linear_binary', obs_params);

% fit
est = tapas_fitModel(sim.y, u, prc_model_config, obs_model_config, optim_config);

% plot
% tapas_hgf_binary_plotTraj(est);


%% Parameter recovery


% because sampleModel doesn't work properly lol
obs_model_config.model = 'tapas_logrt_linear_binary'; 

% set priors
obs_model_config.be2sa = 0; % turn these off
obs_model_config.be3sa = 0;
obs_model_config.be4sa = 0;
obs_model_config.logzesa = 32; % lets be silly
obs_model_config = tapas_align_priors(obs_model_config);


optim_config.nRandInit = 5;


N=200;

% preallocate space
recov.om2.sim = nan(N, 1);
recov.be0.sim = nan(N, 1);
recov.be1.sim = nan(N, 1);
recov.be2.sim = nan(N, 1);
recov.be3.sim = nan(N, 1);
recov.ze.sim = nan(N, 1);
recov.om2.est = nan(N, 1);
recov.be0.est = nan(N, 1);
recov.be1.est = nan(N, 1);
recov.be2.est = nan(N, 1);
recov.be3.est = nan(N, 1);
recov.ze.est = nan(N, 1);
recov.om2.space = 'native';
recov.be0.space = 'native';
recov.be1.space = 'native';
recov.be2.space = 'native';
recov.be3.space = 'native';
recov.ze.space = 'log';
recov.LME = nan(N, 1); % store LME
recov.AIC = nan(N, 1); % store AIC
recov.BIC = nan(N, 1); % store BIC




% Loop
for i = 1:N
    % simulate data with sampleModel
    sim = tapas_sampleModel(u, prc_model_config, obs_model_config);

    % store simulated params
    recov.om2.sim(i) = sim.p_prc.om(2);
    recov.be0.sim(i) = sim.p_obs.be0;
    recov.be1.sim(i) = sim.p_obs.be1;
    recov.be2.sim(i) = sim.p_obs.be2;
    recov.be3.sim(i) = sim.p_obs.be3;
    recov.ze.sim(i) = sim.p_obs.ze;
    
    % recover
    est = tapas_fitModel(...
                    sim.y,...
                    sim.u,...
                    prc_model_config,...
                    obs_model_config,...
                    optim_config);
    
    % store recovered params
    recov.om2.est(i) = est.p_prc.om(2);
    recov.be0.est(i) = est.p_obs.be0;
    recov.be1.est(i) = est.p_obs.be1;
    recov.be2.est(i) = est.p_obs.be2;
    recov.be3.est(i) = est.p_obs.be3;
    recov.ze.est(i) = est.p_obs.ze;
    
    % store fit metrics
    if ~isinf(est.optim.LME)
        recov.LME(i) = est.optim.LME;
        recov.AIC(i) = est.optim.AIC;
        recov.BIC(i) = est.optim.BIC;
    end

end

addpath('..')
recovery_figures(recov)

