function [prc_config, obs_config] = STE_PredictionBias_config

prc_config = prc_PredictionBias_ehgf_binary_pu_tbt_config(); % perceptual model
obs_config = obs_PredictionBias_comb_obs_config(); % response model


prc_config.ommu(2)    = -2;
prc_config.omsa(2)    = 4;

prc_config.rhomu(2)   = 0; % bias towards sad - does work in terms of # responses, but psychometric functions look wrong
prc_config.rhosa(2)   = 4;

prc_config.logalmu    = log(.05); % perceptual uncertainty
prc_config.logalsa    = 2;

prc_config = tapas_align_priors(prc_config);

obs_config.logzemu = log(1);
obs_config.logzesa = 2;

obs_config.beta0mu = 6.5000;
obs_config.beta0sa = 4;

obs_config.beta1mu = 0;
obs_config.beta1sa = 4;

obs_config.beta2mu = 0;
obs_config.beta2sa = 4;

obs_config.beta3mu = 0;
obs_config.beta3sa = 4;

obs_config.beta4mu = 4;
obs_config.beta4sa = 4;

obs_config.logsasa = log(.1);
obs_config.logsasa = 2;

obs_config = tapas_align_priors(obs_config);


end