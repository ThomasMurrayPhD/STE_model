function [c] = obs_ResponseBias_comb_obs_config()
% [c] = m1_comb_obs_config()
%
% Contains the prior configurations of the combined response model M1.
% (Designed to be compatible with the HGF Toolbox as part of TAPAS).
%
% INPUT
%   argin         type  
%
%   OPTIONAL:
%
% OUTPUT    
%   c             struct       Struct containing obs model prior configs
%
% _________________________________________________________________________
% Author: Alex Hess
%
% Copyright (C) 2023 Translational Neuromodeling Unit
%                    Institute for Biomedical Engineering
%                    University of Zurich & ETH Zurich
%
% This file is released under the terms of the GNU General Public Licence
% (GPL), version 3. You can redistribute it and/or modify it under the
% terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any
% later version.
%
% This file is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
% _________________________________________________________________________


%% Config structure
c = struct;

% Is the decision based on predictions or posteriors? Comment as appropriate.
% c.predorpost = 1; % Predictions
c.predorpost = 2; % Posteriors

% Model name
c.model = 'obs_ResponseBias_comb_obs';


%% Sufficient statistics of Gaussian parameter priors
%-----------------------------------
% Model for binary predictions
% zeta0 - inv decision noise
c.logzeta0mu = log(48);
c.logzeta0sa = 1;

% zeta1 - response bias
c.zeta1mu = 0;
c.zeta1sa = 2;

%-----------------------------------
% Model for Reaction Time fit
% beta0 - intercept, mean logrt per participant
c.beta0mu = 6.5;
c.beta0sa = 2;

% beta1 - expected state in "correct space"
c.beta1mu = 0;
c.beta1sa = 2;%2;

% beta2 - sahat 1 in "correct space"
c.beta2mu = 0;
c.beta2sa = 5;%5;

% beta3 - sahat2 
c.beta3mu = 0;
c.beta3sa = 1;%1;

% % beta4 - mu3
c.beta4mu = 0;
c.beta4sa = 0;%1;

% Sigma (noise term)
c.logsamu = -2;
c.logsasa = 2;


%% Gather prior settings in vectors
c.priormus = [
    c.logzeta0mu,...
    c.zeta1mu,...
    c.beta0mu, ...
    c.beta1mu, ...
    c.beta2mu, ...
    c.beta3mu, ...
    c.beta4mu, ...
    c.logsamu ...
         ];

c.priorsas = [
    c.logzeta0sa,...
    c.zeta1sa,...
    c.beta0sa, ...
    c.beta1sa, ...
    c.beta2sa, ...
    c.beta3sa, ...
    c.beta4sa, ...
    c.logsasa ...
    ];

% Model filehandle
c.obs_fun = @obs_ResponseBias_comb_obs;

% Handle to function that transforms perceptual parameters to their native
% space from the space they are estimated in
c.transp_obs_fun = @obs_ResponseBias_comb_obs_transp;

return;
