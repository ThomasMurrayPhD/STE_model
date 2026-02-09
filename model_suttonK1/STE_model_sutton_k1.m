% Script to use the sutton k1 model

% Perceptual model 
%   = [prc3] Sutton K1
% Response model 
%   = [obs4] combined unitsq_sgm and logRT (sutton_k1 specific)


%%
clear;
close all;
addpath('..');

%%

% example data (to get contingencies etc)
sub_data = readtable('..\STE_data\10369536_A_Threat.csv');

% Contingency space
cue = sub_data.Cue_idx;
cue(cue==0) = -1; % balanced contrast coding for the "happy bias"
outcome = sub_data.Outcome_p_sad/100;
u_al = 1 - (0.5*(1+cue)-cue.*outcome);
state = double(sub_data.Cue_idx == sub_data.Outcome_idx);

sub_data.u_al = u_al;
sub_data.state=state;
sub_data.p_sad=outcome;


% Responses
sub_data.logRT = log(sub_data.Response_RT);
sub_data.resp_state = double(sub_data.Cue_idx == sub_data.Response_idx);

% model input
% u = [sub_data.u_al, cue];
% y = [sub_data.resp_state, sub_data.logRT];%, sub_data.Confidence_idx];

u = sub_data.state;
y = sub_data.resp_state;


%%
prc_model_config = prc_sutton_k1_binary_config();
prc_model_config.logitvhat_1sa = 0; % fix vhat1 to .5
prc_model_config = tapas_align_priors(prc_model_config);

bopars = tapas_fitModel([],...
                         u,...
                         prc_model_config,...
                         'tapas_bayes_optimal_binary_config',...
                         'tapas_quasinewton_optim_config');
%%
% 
% % Ok so sutton_k1 doesn't work with posteriors (only predictions)...
% % obs_model_config.predorpost = 1;
% % 
% % Maybe I could tru editing the infStates in tapas_sutton_k1 - at the
% % moment it only outputs vhat, but not v
% % 
% % ^This is what I have done.
% 
% 
% sim = tapas_simModel(u,...
%                      'prc_sutton_k1_binary',...
%                      bopars.p_prc.p,...
%                      'obs_suttonK1_unitsq_sgm',...
%                      5,...
%                      123456789);
% 
% % prc_sutton_k1_binary_plotTraj(sim);
% 
% % Fit
% est = tapas_fitModel(...
%     sim.y,...
%     sim.u,...
%     prc_model_config,...
%     obs_model_config,...
%     optim_config);
% 
% 
% % Cool... Ok lets try out the logRT model
% 
% %% LogRT 
% prc_model_config = prc_sutton_k1_binary_config();
% obs_model_config = obs_suttonK1_logrt_linear_binary_config();
% optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
% 
% prc_model_config.logitvhat_1sa = 0; % fix vhat1 to .5
% prc_model_config = tapas_align_priors(prc_model_config);
% 
% obs_params = obs_suttonK1_logrt_linear_binary_transp([], obs_model_config.priormus);
% 
% sim = tapas_simModel(u,...
%                      'prc_sutton_k1_binary',...
%                      bopars.p_prc.p,...
%                      'obs_suttonK1_logrt_linear_binary',...
%                      obs_params,...
%                      123456789);
% 
% est = tapas_fitModel(...
%     sim.y,...
%     sim.u,...
%     prc_model_config,...
%     obs_model_config,...
%     optim_config);
% % Seems to work....


%% Combined
% Ok let's try the combined response model
prc_model_config = prc_sutton_k1_binary_config();

prc_model_config.logmumu = log(3); % set mu

prc_model_config.logitvhat_1sa = tapas_logit(0.5, 1); % fix vhat1 to .5
prc_model_config.logitvhat_1sa = 0; % fix vhat1 to .5

prc_model_config.logRhatmu = log(.5); % set Rhat
prc_model_config.logh_1mu = log(.005);

prc_model_config = tapas_align_priors(prc_model_config);



obs_model_config = obs_suttonK1_comb_obs_config();
obs_model_config.beta1mu = -1;
obs_model_config.beta2mu = 1;
obs_model_config = tapas_align_priors(obs_model_config);

optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 10; % Annoying but fits better

prc_params = prc_sutton_k1_binary_transp([], prc_model_config.priormus);
obs_params = obs_suttonK1_comb_obs_transp([], obs_model_config.priormus);

sim = tapas_simModel(u,...
                     'prc_sutton_k1_binary',...
                     prc_params,...
                     'obs_suttonK1_comb_obs',...
                     obs_params,...
                     123456789);

% plot trajectory
prc_sutton_k1_binary_plotTraj(sim);


est = tapas_fitModel(...
    sim.y,...
    sim.u,...
    prc_model_config,...
    obs_model_config,...
    optim_config);




%% Full parameter recovery

% set perceptual priors
prc_model_config = prc_sutton_k1_binary_config();
prc_model_config.logmumu = log(3); % set mu prior mean
prc_model_config.logmusa = 16; % set mu prior mean
prc_model_config.logRhatmu = log(.5); % set prior mean Rhat
prc_model_config.logRhatsa = 16; % set prior mean Rhat
prc_model_config.logitvhat_1sa = tapas_logit(0.5, 1); % set vhat1 to .5
prc_model_config.logitvhat_1sa = 0; % fix
prc_model_config.logh_1mu = log(.005);
prc_model_config.logh_1sa = 4;
prc_model_config = tapas_align_priors(prc_model_config);

% set obs model priors
obs_model_config = obs_suttonK1_comb_obs_config();
obs_model_config.logzemu = log(48);
obs_model_config.logzesa = 8;
obs_model_config.beta1mu = -1;
obs_model_config.beta1sa = 8;
obs_model_config.beta2mu = 1;
obs_model_config.beta2sa = 8;
obs_model_config.logsamu = log(3);
obs_model_config.logsasa = 8;
obs_model_config = tapas_align_priors(obs_model_config);

% optimisation algorithm
optim_config     = tapas_quasinewton_optim_config();
optim_config.nRandInit = 10; % Annoying but fits better

% number of iterations
N=200;

% Parameters to recover
prc_param_names = {'mu', 'Rhat', 'h_1'};
prc_param_idx   = [1, 2, 4];
prc_param_space = {'log', 'log', 'log'};
obs_param_names = {'ze', 'beta0', 'beta1', 'beta2', 'sa'};
obs_param_idx   = [1, 2, 3, 4, 5];
obs_param_space = {'log', 'native', 'native', 'native', 'log'};


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
save('model_suttonK1_recovery.mat', 'recov');
recovery_figures(recov);


% Check parameter values that result in LME=-Inf








