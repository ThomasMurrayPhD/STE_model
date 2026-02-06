function [logp, yhat, res, logp_split] = obs_suttonK1_comb_obs(r, infStates, ptrans)
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

% Zeta parameter (for binary) - check exp
ze = exp(ptrans(1));

% Parameters for logRT - check transformations
be0  = ptrans(2);
be1  = ptrans(3);
be2  = ptrans(4);
sa   = exp(ptrans(5));


%% binary part of the response model
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
x = infStates(:,1,pop);
x(r.irr) = [];
y = r.y(:,1);
y(r.irr) = [];

% Avoid any numerical problems when taking logarithms close to 1
logx = log(x);
log1pxm1 = log1p(x-1);
logx(1-x<1e-4) = log1pxm1(1-x<1e-4);
log1mx = log(1-x);
log1pmx = log1p(-x);
log1mx(x<1e-4) = log1pmx(x<1e-4); 

% Calculate log-probabilities for non-irregular trials
reg = ~ismember(1:n,r.irr);
logp_binary(reg) = y.*ze.*(logx -log1mx) +ze.*log1mx -log((1-x).^ze +x.^ze);
yhat_binary(reg) = x;
res_binary(reg) = (y-x)./sqrt(x.*(1-x));


%% continuous part of the response model
% Transform parameters to their native space


% Initialize returned log-probabilities, predictions,
% and residuals as NaNs so that NaN is returned for all
% irregualar trials
n = size(infStates,1);
logp_reactionTime = NaN(n,1);
yhat_reactionTime = NaN(n,1);
res_reactionTime  = NaN(n,1);

% Weed irregular trials out from responses and inputs
y = r.y(:,1);
y(r.irr) = [];

u = r.u(:,1);
u(r.irr) = [];

% Extract trajectories of interest from infStates
da = r.traj.da; % prediction error
be = r.traj.be; % beta - unconstrained learning rate (log gain)
al = r.traj.al; % alpha
h = r.traj.h; % h (not sure)
v = r.traj.v; % posterior
vhat = r.traj.vhat; % prediction



% Calculate predicted log-reaction time
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
logrt = be0 +be1.*al +be2.*abs(da);

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
