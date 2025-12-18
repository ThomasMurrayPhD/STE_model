% Script to use the sutton k1 model

% Perceptual model 
%   = [prc2] Sutton K1
% Response model 
%   = [obs1] combined unitsq_sgm and logRT


%%
clear;
close all;

%%

% example data (to get contingencies etc)
sub_data = readtable('STE_data\10369536_A_Threat.csv');

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
prc_model_config = tapas_sutton_k1_binary_config();
obs_model_config = tapas_unitsq_sgm_config();

bopars = tapas_fitModel([],...
                         u,...
                         'tapas_sutton_k1_binary_config',...
                         'tapas_bayes_optimal_binary_config',...
                         'tapas_quasinewton_optim_config');
%%

% Ok so sutton_k1 doesn't work with posteriors (only predictions)...
% obs_model_config.predorpost = 1;

sim = tapas_simModel(u,...
                     'tapas_sutton_k1_binary',...
                     bopars.p_prc.p,...
                     'tapas_unitsq_sgm',...
                     5,...
                     123456789);


tapas_sutton_k1_binary_plotTraj(sim)


