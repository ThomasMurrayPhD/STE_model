% Function to run parameter recovery for PerceptualBias model


% To do
%   - check sampleModel is using correct model - need to specify correct
%   model name in field in obs configs
%   - Maybe edit sampleModel to spit out the parameters that it sampled,
%   regardless of whether it crashes the model or not.
%   - Create config files




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
y = [sub_data.resp_state, sub_data.logRT];%, sub_data.Confidence_idx];


%% Get configuration structures
[prc_config, obs_config] = STE_PerceptualBias_config;
optim_config     = tapas_quasinewton_optim_config(); % optimisation algorithm
optim_config.nRandInit = 5;


%% Run parameter recovery

% run recovery
N=200;

% preallocate space
recov.rho.sim = nan(N, 1);
recov.rho.est = nan(N, 1);
recov.rho.space = 'native';
recov.om2.sim = nan(N, 1);
recov.om2.est = nan(N, 1);
recov.om2.space = 'native';
recov.al.sim = nan(N, 1);
recov.al.est = nan(N, 1);
recov.al.space = 'log';
recov.ze.sim = nan(N, 1);
recov.ze.est = nan(N, 1);
recov.ze.space = 'log';
recov.beta0.sim = nan(N, 1);
recov.beta0.est = nan(N, 1);
recov.beta0.space = 'native';
recov.beta1.sim = nan(N, 1);
recov.beta1.est = nan(N, 1);
recov.beta1.space = 'native';
recov.beta2.sim = nan(N, 1);
recov.beta2.est = nan(N, 1);
recov.beta2.space = 'native';
recov.beta3.sim = nan(N, 1);
recov.beta3.est = nan(N, 1);
recov.beta3.space = 'native';
recov.beta4.sim = nan(N, 1);
recov.beta4.est = nan(N, 1);
recov.beta4.space = 'native';
recov.sa.sim = nan(N, 1);
recov.sa.est = nan(N, 1);
recov.sa.space = 'log';

recov.LME = nan(N, 1); % store LME
recov.AIC = nan(N, 1); % store AIC
recov.BIC = nan(N, 1); % store BIC

% Loop
for i = 1:N
    
    try
        % simulate data with sampleModel
        sim = tapas_sampleModel(u, prc_config, obs_config);

        % store simulated params
        recov.rho.sim(i)    = sim.p_prc.rho;
        recov.om2.sim(i)    = sim.p_prc.om(2);
        recov.al.sim(i)     = sim.p_prc.al;
        recov.ze.sim(i)     = sim.p_obs.ze;
        recov.beta0.sim(i)  = sim.p_obs.beta0;
        recov.beta1.sim(i)  = sim.p_obs.beta1;
        recov.beta2.sim(i)  = sim.p_obs.beta2;
        recov.beta3.sim(i)  = sim.p_obs.beta3;
        recov.beta4.sim(i)  = sim.p_obs.beta4;
        recov.sa.sim(i)     = sim.p_obs.sa;

        % recover
        est = tapas_fitModel(...
                    sim.y,...
                    sim.u,...
                    prc_config,...
                    obs_config,...
                    optim_config);
    
        % store recovered params
        recov.rho.est(i)    = est.p_prc.rho;
        recov.om2.est(i)    = est.p_prc.om(2);
        recov.al.est(i)     = est.p_prc.al;
        recov.ze.est(i)     = est.p_obs.ze;
        recov.beta0.est(i)  = est.p_obs.beta0;
        recov.beta1.est(i)  = est.p_obs.beta1;
        recov.beta2.est(i)  = est.p_obs.beta2;
        recov.beta3.est(i)  = est.p_obs.beta3;
        recov.beta4.est(i)  = est.p_obs.beta4;
        recov.sa.est(i)     = est.p_obs.sa;

        % store fit metrics
        if ~isinf(est.optim.LME)
            recov.LME(i) = est.optim.LME;
            recov.AIC(i) = est.optim.AIC;
            recov.BIC(i) = est.optim.BIC;
        end

        % Store full est
        recov.est_full{i} = est; % This might crash the figures plot
    catch
    end

end

save('PerceptualBias_recovery.mat', 'recov');
recovery_figures(recov);


