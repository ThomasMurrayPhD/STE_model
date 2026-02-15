% Function to run parameter recovery for PerceptualBias model

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



%% Run parameter recovery

% % run recovery
% N=200;
% 
% % Parameters to recover
% prc_param_names = {'om2', 'rho2', 'al'};
% prc_param_idx   = [13, 8, 15];
% prc_param_space = {'native', 'native', 'log'};
% obs_param_names = {'ze', 'beta0', 'beta1', 'beta2', 'beta3', 'beta4', 'sa'};
% obs_param_idx   = [1, 2, 3, 4, 5, 6, 7];
% obs_param_space = {'log', 'native', 'native', 'native', 'native', 'native', 'log'};
% 
% % preallocate sim and est parameters
% all_params = [prc_param_names, obs_param_names];
% for iP = 1:numel(all_params)
%     recov.(all_params{iP}).sim = nan(N, 1);
%     recov.(all_params{iP}).est = nan(N, 1);
% end
% recov.LME = nan(N, 1); % store LME
% recov.AIC = nan(N, 1); % store AIC
% recov.BIC = nan(N, 1); % store BIC
% recov.est = cell(N,1); % store full est structure
% 
% % Main loop
% for i = 1:N
% 
%     % Sample from model params and simulate
%     sim = tapas_sampleModel(u, prc_config, obs_config);
% 
%     % Store simulated prc params
%     for iP = 1:numel(prc_param_names)
%         param_name=prc_param_names{iP};
%         recov.(param_name).sim(i) = sim.p_prc.p(prc_param_idx(iP));
%         recov.(param_name).space = prc_param_space{iP};
%     end
% 
%     % store simulated obs params
%     for iP = 1:numel(obs_param_names)
%         param_name=obs_param_names{iP};
%         recov.(param_name).sim(i) = sim.p_obs.p(obs_param_idx(iP));
%         recov.(param_name).space = obs_param_space{iP};
%     end
% 
% 
%     % estimate parameters from simulated responses
%     if ~any(isnan(sim.y)) % only if no missing trials
%         try
%             est = tapas_fitModel(...
%                 sim.y,...
%                 sim.u,...
%                 prc_config,...
%                 obs_config,...
%                 optim_config);
% 
%             % get estimated parameters
%             if ~isinf(est.optim.LME)
%                 % store fit metrics
%                 recov.LME(i) = est.optim.LME;
%                 recov.AIC(i) = est.optim.AIC;
%                 recov.BIC(i) = est.optim.BIC;
% 
%                 % Store estimated prc params
%                 for iP = 1:numel(prc_param_names)
%                     param_name=prc_param_names{iP};
%                     recov.(param_name).est(i) = est.p_prc.p(prc_param_idx(iP));
%                 end
% 
%                 % store simulated obs params
%                 for iP = 1:numel(obs_param_names)
%                     param_name=obs_param_names{iP};
%                     recov.(param_name).est(i) = est.p_obs.p(obs_param_idx(iP));
%                 end
% 
%             end
%             recov.est{i} = est;
%         catch
%             fprintf('\nCannot fit')
%         end
%     end
% end




%% Recovery 1
addpath('..')
% Hmm so this isn't actually very good. Lots of bad fits in beta0 and
% log(ze)
% save('PredictionBias_recovery2.mat', 'recov');
% recovery_figures(recov);
% 
% all_params = fieldnames(recov);
% all_params(ismember(all_params, {'LME', 'AIC', 'BIC', 'est'})) = [];
% 
% bad_fits = recov.beta0.est==6.5;
% 
% for iP = 1:numel(all_params)
%     sim = recov.(all_params{iP}).sim;
%     est = recov.(all_params{iP}).est;
%     valid = ~isnan(est);
%     sim = sim(valid);
%     est = est(valid);
% 
%     bad_fits_param = bad_fits(valid);
% 
%     figname = all_params{iP};
%     if strcmp(recov.(all_params{iP}).space, 'log')
%         sim = log(sim);
%         est = log(est);
%         figname = ['log(', figname, ')'];
%     elseif strcmp(recov.(all_params{iP}).space, 'logit')
%         sim = tapas_logit(sim, 1);
%         est = tapas_logit(est, 1);
%         figname = ['logit(', figname, ')'];
%     end
% 
%     figure('name', figname);  
%     hold on;
%     scatter(sim, est, 40);
%     scatter(sim(bad_fits_param), est(bad_fits_param), 40, 'red');
%     h = refline(1,0); 
%     h.Color = [0.8500    0.3250    0.0980];
%     xlabel('Simulated');
%     ylabel('Estimated');
%     corr_string = sprintf("%s\nPearson's r = %.3f; Spearman's \\rho = %.3f; Kendall's \\tau = %.3f", ...
%         figname,...
%         corr(sim, est), ...
%         corr(sim, est, 'type', 'Spearman'), ...
%         corr(sim, est, 'type', 'Kendall'));
%     title(corr_string)
% end
% 
% 



%% Change priors...
% If this works, update STE_PredictionBias_config;
% 
% If it doesn't - maybe use stim_noise (not pu) in beta4 in this model?
obs_config.beta1sa = 2;
obs_config = tapas_align_priors(obs_config);



%% Run parameter recovery again...


% run recovery
N=200;

% Parameters to recover
prc_param_names = {'om2', 'rho2', 'al'};
prc_param_idx   = [13, 8, 15];
prc_param_space = {'native', 'native', 'log'};
obs_param_names = {'ze', 'beta0', 'beta1', 'beta2', 'beta3', 'beta4', 'sa'};
obs_param_idx   = [1, 2, 3, 4, 5, 6, 7];
obs_param_space = {'log', 'native', 'native', 'native', 'native', 'native', 'log'};

% preallocate sim and est parameters
all_params = [prc_param_names, obs_param_names];
for iP = 1:numel(all_params)
    recov.(all_params{iP}).sim = nan(N, 1);
    recov.(all_params{iP}).est = nan(N, 1);
end
recov.LME = nan(N, 1); % store LME
recov.AIC = nan(N, 1); % store AIC
recov.BIC = nan(N, 1); % store BIC
recov.est = cell(N,1); % store full est structure

% Main loop
for i = 1:N
    
    % Sample from model params and simulate
    sim = tapas_sampleModel(u, prc_config, obs_config);
    
    % Store simulated prc params
    for iP = 1:numel(prc_param_names)
        param_name=prc_param_names{iP};
        recov.(param_name).sim(i) = sim.p_prc.p(prc_param_idx(iP));
        recov.(param_name).space = prc_param_space{iP};
    end
    
    % store simulated obs params
    for iP = 1:numel(obs_param_names)
        param_name=obs_param_names{iP};
        recov.(param_name).sim(i) = sim.p_obs.p(obs_param_idx(iP));
        recov.(param_name).space = obs_param_space{iP};
    end


    % estimate parameters from simulated responses
    if ~any(isnan(sim.y)) % only if no missing trials
        try
            est = tapas_fitModel(...
                sim.y,...
                sim.u,...
                prc_config,...
                obs_config,...
                optim_config);

            % get estimated parameters
            if ~isinf(est.optim.LME)
                % store fit metrics
                recov.LME(i) = est.optim.LME;
                recov.AIC(i) = est.optim.AIC;
                recov.BIC(i) = est.optim.BIC;
                
                % Store estimated prc params
                for iP = 1:numel(prc_param_names)
                    param_name=prc_param_names{iP};
                    recov.(param_name).est(i) = est.p_prc.p(prc_param_idx(iP));
                end
                
                % store simulated obs params
                for iP = 1:numel(obs_param_names)
                    param_name=obs_param_names{iP};
                    recov.(param_name).est(i) = est.p_obs.p(obs_param_idx(iP));
                end

            end
            recov.est{i} = est;
        catch
            fprintf('\nCannot fit')
        end
    end
end
save('PredictionBias_recovery3.mat', 'recov');
recovery_figures(recov);








