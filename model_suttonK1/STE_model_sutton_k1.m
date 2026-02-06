% Script to use the sutton k1 model

% Perceptual model 
%   = [prc3] Sutton K1
% Response model 
%   = [obs4] combined unitsq_sgm and logRT (sutton_k1 specific)


%%
% clear;
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
obs_model_config = obs_suttonK1_unitsq_sgm_config();
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm

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

%% LogRT 
% prc_model_config = prc_sutton_k1_binary_config();
% obs_model_config = obs_suttonK1_logrt_linear_binary_config();
% optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
% 
% prc_model_config.logitvhat_1sa = 0; % fix vhat1 to .5
% prc_model_config = tapas_align_priors(prc_model_config);
% 
% 
% sim = tapas_simModel(u,...
%                      'prc_sutton_k1_binary',...
%                      bopars.p_prc.p,...
%                      'obs_suttonK1_logrt_linear_binary',...
%                      obs_model_config.priormus,...
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
prc_model_config.logitvhat_1sa = 0; % fix vhat1 to .5
prc_model_config = tapas_align_priors(prc_model_config);

obs_model_config = obs_suttonK1_comb_obs_config();
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 15; % Annoying but fits better

sim = tapas_simModel(u,...
                     'prc_sutton_k1_binary',...
                     bopars.p_prc.p,...
                     'obs_suttonK1_comb_obs',...
                     obs_model_config.priormus,...
                     123456789);

est = tapas_fitModel(...
    sim.y,...
    sim.u,...
    prc_model_config,...
    obs_model_config,...
    optim_config);




