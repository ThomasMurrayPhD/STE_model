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
[prc_config, obs_config] = STE_RW_config;

bopars = tapas_fitModel([],...
                         u,...
                         'prc_RW_binary_config',...
                         'tapas_bayes_optimal_binary_config',...
                         'tapas_quasinewton_optim_config');
%%


sim = tapas_simModel(u,...
                     'prc_RW_binary',...
                     bopars.p_prc.p,...
                     'obs_RW_unitsq_sgm',...
                     5,...
                     123456789);

prc_RW_binary_plotTraj(sim);


