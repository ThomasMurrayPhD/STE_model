function [logp, yhat, res, logp_split] = obs_PerceptualBias_comb_obs(r, infStates, ptrans)
% [logp, yhat, res, logp_split] = m1_comb_obs(r, infStates, ptrans)
%
% Calculates the combined log-probability of binary and continuous
% behavioural responses.
%
% INPUT
%   r             struct      Struct obtained from tapas_fitModel.m fct
%   infStates     tensor      Tensor containing inferred states from the
%                             perceptual model    
%   ptrans        vector      1xP vector with free param values (est space)
%
%   OPTIONAL:
%
% OUTPUT    
%   logp          vector       1xN vector containing trialwise log
%                              probabilities of responses
%   yhat          matrix       Nx2 matrix containing noise-free predictions
%   res           matrix       Nx2 matrix containing responses (bin + cont)
%   logp_split    matrix       Nx2 matrix containing trialwise logll values
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


%% Separate parameters
% Transform parameters to their native space

% Parameters for binary
ze = exp(ptrans(1));

% Parameters for continuous
be0  = ptrans(2);
be1  = ptrans(3);
be2  = ptrans(4);
be3  = ptrans(5);
be4  = ptrans(6);
sa   = exp(ptrans(7));

%% binary part of the response model
% compute log likelihood (binary responses)

% Predictions or posteriors?
pop = 1; % Default: predictions
if r.c_obs.predorpost == 2
    pop = 3; % Alternative: posteriors
end

% Initialize returned log-probabilities as NaNs so that NaN is
% returned for all irregualar trials
n = size(infStates,1);
logp_binary = NaN(n,1);
yhat_binary = NaN(n,1);
res_binary  = NaN(n,1);

% Weed irregular trials out from inferred states and responses
x_state = infStates(:,1,pop);

x_state(r.irr) = [];
y = r.y(:,1);
y(r.irr) = [];

% Avoid any numerical problems when taking logarithms close to 1
logx = log(x_state);
log1pxm1 = log1p(x_state-1);
logx(1-x_state<1e-4) = log1pxm1(1-x_state<1e-4);
log1mx = log(1-x_state);
log1pmx = log1p(-x_state);
log1mx(x_state<1e-4) = log1pmx(x_state<1e-4); 


% Calculate log-probabilities for non-irregular trials
reg = ~ismember(1:n,r.irr);
logp_binary(reg) = y.*ze.*(logx -log1mx) +ze.*log1mx -log((1-x_state).^ze +x_state.^ze);
yhat_binary(reg) = x_state;
res_binary(reg) = (y-x_state)./sqrt(x_state.*(1-x_state));




%% continuous part of the response model
% Compute the log likelihood (logRTs)


% Initialize returned log-probabilities, predictions,
% and residuals as NaNs so that NaN is returned for all
% irregualar trials
n = size(infStates,1);
logp_reactionTime = NaN(n,1);
yhat_reactionTime = NaN(n,1);
res_reactionTime  = NaN(n,1);

% Weed irregular trials out from responses and inputs
y = r.y(:,2);
y(r.irr) = [];


u_al = r.u(:,1);
u = u_al>0.5;
stim_noise = 0.5-abs(u_al-.5); % [0,1]->0, [.2,.8]->.2, [.4,.6]->.4


% Extract trajectories of interest from infStates
mu1hat = infStates(:,1,1);
sa1hat = infStates(:,1,2);
mu2    = infStates(:,2,3);
sa2    = infStates(:,2,4);
% mu3    = infStates(:,3,3);


% Surprise
% ~~~~~~~~
poo = mu1hat.^u.*(1-mu1hat).^(1-u); % probability of observed outcome
surp = -log2(poo);
% surp_shifted = [1; surp(1:(length(surp)-1))];

% Bernoulli variance (aka irreducible uncertainty, risk) 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bernv = sa1hat;

% Inferential variance (aka informational or estimation uncertainty, ambiguity)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inferv = tapas_sgm(mu2, 1).*(1 -tapas_sgm(mu2, 1)).*sa2; % transform down to 1st level

% Phasic volatility (aka environmental or unexpected uncertainty)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% pv = tapas_sgm(mu2, 1).*(1-tapas_sgm(mu2, 1)).*exp(mu3); % transform down to 1st level

mu1 = infStates(:,1,3);
pu = 0.5-abs(mu1-0.5);

% Calculate predicted log-reaction time
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% logrt = be0 +be1.*surp +be2.*bernv +be3.*inferv +be4.*pv;
% logrt = be0 +be1.*surp +be2.*bernv +be3.*inferv +be4.*stim_noise;
logrt = be0 +be1.*surp +be2.*bernv +be3.*inferv +be4.*pu; %%% THIS SEEMS
% TO WORK... MAX RT SHIFTS WITH PSE


% Remove missed/"irregular" trials
logrt(r.irr) = [];


% Calculate log-probabilities for non-irregular trials
% Note: 8*atan(1) == 2*pi (this is used to guard against
% errors resulting from having used pi as a variable).
reg = ~ismember(1:n,r.irr);
logp_reactionTime(reg) = -1/2.*log(8*atan(1).*sa) -(y-logrt).^2./(2.*sa);
yhat_reactionTime(reg) = logrt;
res_reactionTime(reg) = y-logrt;



%% get combined log likelihood of two response data modalities
logp = logp_binary + logp_reactionTime;
logp_split = [logp_binary logp_reactionTime];

%% return predictions and responses for each response data modality
yhat = [yhat_binary yhat_reactionTime];
res = [res_binary res_reactionTime];

end
