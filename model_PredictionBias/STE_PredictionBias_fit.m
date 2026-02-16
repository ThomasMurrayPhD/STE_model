% Function to run parameter recovery for PredictionBias model

clear; close all; clc


%% Get inputs
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
y = [sub_data.resp_state, sub_data.logRT];


%% Get configuration structures
[prc_config, obs_config] = STE_PredictionBias_config;
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;


%% Fit data
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

    % subject responses
    sub_data.logRT = log(sub_data.Response_RT);
    sub_data.resp_state = double(sub_data.Cue_idx == sub_data.Response_idx);
    y = [sub_data.resp_state, sub_data.logRT];

    % remove missing
    missed = isnan(sub_data.Response_idx);
    y(missed, :) = NaN;

    % fit model
    try
        model_fits(i).est = tapas_fitModel(...
            y,...
            u,...
            prc_config,...
            obs_config,...
            optim_config);
    catch
    end
    completion_times(i) = toc;
end

save('model_PredictionBias_fit2.mat', 'model_fits');
