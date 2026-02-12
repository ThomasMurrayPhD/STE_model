%% MODEL 3
% Perceptual model 
%   = [prc1] 2 level eHGF Nace remix (fixed rho=0)
% Response model 
%   = [obs2] combined unitsq_sgm and logRT (bias on response)


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
u = [sub_data.u_al, cue];
y = [sub_data.resp_state, sub_data.logRT];%, sub_data.Confidence_idx];


%% Get configuration structures

[prc_config, obs_config] = STE_ResponseBias_config;
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;


%% simulate responses

r_temp = [];
r_temp.c_prc.n_levels = 3;
prc_params = prc_ResponseBias_ehgf_binary_pu_tbt_transp(r_temp, prc_config.priormus);

obs_params = obs_config.priormus;
obs_params(1) = exp(obs_params(1));
obs_params(8) = exp(obs_params(8));

sim_3_bias = tapas_simModel(u,...
    prc_config.model,...
    prc_params,...
    obs_config.model,...
    obs_params,...
    123456789);


sim_sad = (sub_data.Cue_idx == 1 & sim_3_bias.y(:,1) == 1) + (sub_data.Cue_idx == 0 & sim_3_bias.y(:,1) == 0);
N_sad = sum(sim_sad);

visualise_psychometric(u, sub_data, 'prc_ResponseBias_ehgf_binary_pu_tbt', prc_params, 'obs_ResponseBias_comb_obs', obs_params, 20)


%% Plot trajectory
prc_ResponseBias_ehgf_binary_tbt_plotTraj(sim_3_bias);


%% recover parameters
est = tapas_fitModel(...
    sim_3_bias.y,...
    sim_3_bias.u,...
    prc_model_config,...
    obs_model_config,...
    optim_config);





