function [prc_config, obs_config] = STE_ResponseBias_config

prc_config = prc_ResponseBias_ehgf_binary_pu_tbt_config(); % perceptual model
obs_config = obs_ResponseBias_comb_obs_config(); % response model

prc_config.ommu(2)    = -2;
prc_config.omsa(2)    = 4;

prc_config.rhomu(2)   = 0; % bias towards sad
prc_config.rhosa(2)   = 0; %% IMPORTANT TO FIX IN THIS SETUP

prc_config.logalmu    = log(.1); % perceptual uncertainty
prc_config.logalsa    = 2;

prc_config = tapas_align_priors(prc_config);

obs_config.logzeta0mu = log(1);
obs_config.logzeta0sa = 2;

obs_config.zeta1mu = 3; % Response bias
obs_config.zeta1sa = 2;

obs_config.beta0mu = 6.5000;
obs_config.beta0sa = 4;

obs_config.beta1mu = 0;
obs_config.beta1sa = 4;

obs_config.beta2mu = 0;
obs_config.beta2sa = 4;

obs_config.beta3mu = 0;
obs_config.beta3sa = 4;

obs_config.beta4mu = 2;
obs_config.beta4sa = 4;

obs_config.logsasa = log(.1);
obs_config.logsasa = 2;

obs_config = tapas_align_priors(obs_config);



end