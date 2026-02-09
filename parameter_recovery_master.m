function recov = parameter_recovery_master(u,...
    prc_model_config,...
    obs_model_config,...
    optim_config,...
    N,...
    prc_params,...
    prc_param_idx,...
    prc_param_space,...
    obs_params,...
    obs_param_idx,...
    obs_param_space)

% master function to perform parameter recovery
% 
% input:
%   u                   - model input
%   prc_model_config
%   obs_model_config
%   optim_config
%   N                   - N simulations
%   prc_params          - cell with names of perceptual model params
%   prc_param_idx       - array with idxs of params in prc_params
%   obs_params          - cell with names of observation model params
%   obs_param_idx       - array with idxs of params in obs_params




% preallocate sim and est parameters
all_params = [prc_params, obs_params];
for iP = 1:numel(all_params)
    recov.(all_params{iP}).sim = nan(N, 1);
    recov.(all_params{iP}).est = nan(N, 1);
end
recov.LME = nan(N, 1); % store LME
recov.AIC = nan(N, 1); % store AIC
recov.BIC = nan(N, 1); % store BIC

% Just to be fancy
completion_times = zeros(N, 1);

% Main loop
for i = 1:N
    tic;
    fprintf('\n--------------------------------\n')
    fprintf('\nParameter recovery iteration %i\n', i);
    if i>1
        avg_iter_time = mean(completion_times(1:i));
        fprintf('\tAverage iteration time = %1.2fs', avg_iter_time)
        estimated_total_time = avg_iter_time * ((N-i) + 1);
        fprintf('\n\tEstimated completion time = %im, %1.2fs\n\n', floor(estimated_total_time/60), rem(estimated_total_time,60));
    end

    sim = tapas_sampleModel(u, prc_model_config, obs_model_config);
    
    % Store simulated prc params
    for iP = 1:numel(prc_params)
        param_name=prc_params{iP};
        recov.(param_name).sim(i) = sim.p_prc.p(prc_param_idx(iP));
        recov.(param_name).space = prc_param_space{iP};
    end
    
    % store simulated obs params
    for iP = 1:numel(obs_params)
        param_name=obs_params{iP};
        recov.(param_name).sim(i) = sim.p_obs.p(obs_param_idx(iP));
        recov.(param_name).space = obs_param_space{iP};
    end


    % estimate parameters from simulated responses
    if ~any(isnan(sim.y)) % only if no missing trials
        try
            est = tapas_fitModel(...
                sim.y,...
                sim.u,...
                prc_model_config,...
                obs_model_config,...
                optim_config);

            % get estimated parameters
            if ~isinf(est.optim.LME)
                % store fit metrics
                recov.LME(i) = est.optim.LME;
                recov.AIC(i) = est.optim.AIC;
                recov.BIC(i) = est.optim.BIC;
                
                % Store estimated prc params
                for iP = 1:numel(prc_params)
                    param_name=prc_params{iP};
                    recov.(param_name).est(i) = est.p_prc.p(prc_param_idx(iP));
                end
                
                % store simulated obs params
                for iP = 1:numel(obs_params)
                    param_name=obs_params{iP};
                    recov.(param_name).est(i) = est.p_obs.p(obs_param_idx(iP));
                end

            end
        catch
            fprintf('\nCannot fit')
        end
        completion_times(i) = toc;
    end
end



end