function [y, yhat] = obs_suttonK1_comb_obs_sim(r, infStates, p)
% [y, yhat] = m1_comb_obs_sim(r, infStates, p)
%
% Simulates responses for binary and continuous data modality.
% (Designed to be compatible with the HGF Toolbox as part of TAPAS).
%
% INPUT
%   r             struct      Struct obtained from tapas_simModel.m fct
%   infStates     tensor      Tensor containing inferred states from the
%                             perceptual model    
%   p             vector      1xP vector with free param values (nat space)
%
%   OPTIONAL:
%
% OUTPUT    
%   y             matrix       Nx2 matrix with simulated responses
%   yhat          matrix       Nx2 matrix containing noise-free predictions
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

% Zeta parameter (binary response)
ze = p(1);

% logRT parameters
be0  = p(2);
be1  = p(3);
be2  = p(4);
sa   = p(5);


%% Run sim for binary predictions
% Predictions or posteriors?
pop = 1; % Default: predictions
if r.c_obs.predorpost == 2
    pop = 3; % Alternative: posteriors
end


% Assumed structure of infStates:
% dim 1: time (ie, input sequence number)
% dim 2: HGF level
% dim 3: 1: muhat, 2: sahat, 3: mu, 4: sa

% Belief trajectories at 1st level
states = squeeze(infStates(:,1,pop));

% Apply the unit-square sigmoid to the inferred states
yhat_pred = states.^ze./(states.^ze+(1-states).^ze);

% Initialize random number generator
if isnan(r.c_sim.seed)
    rng('shuffle');
else
    rng(r.c_sim.seed);
end

% Simulate
y_binary = binornd(1, yhat_pred);



%% Run sim for continuous data modality (logRTs)


% Number of trials
n = size(infStates,1);

% Inputs
u = r.u(:,1);

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

% Initialize random number generator
if isnan(r.c_sim.seed)
    rng('shuffle');
else
    rng(r.c_sim.seed);
end

% Simulate
y_reactionTime = logrt+sqrt(sa)*randn(n, 1);
yhat_reactionTime = logrt;


%% save values for both response data modalities
y = [y_binary y_reactionTime];
yhat = [yhat_pred yhat_reactionTime];

end

