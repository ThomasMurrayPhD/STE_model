function [y, yhat] = obs_PerceptualBias_comb_obs_sim(r, infStates, p)
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


% Parameters for binary
ze = p(1);

% Parameters for logRT
be0  = p(2);
be1  = p(3);
be2  = p(4);
be3  = p(5);
be4  = p(6);
sa   = p(7);


%% Run sim for binary predictions

% Predictions or posteriors?
pop = 1; % Default: predictions
if r.c_obs.predorpost == 2
    pop = 3; % Alternative: posteriors
end
x_state = infStates(:,1,pop);


% Assumed structure of infStates:
% dim 1: time (ie, input sequence number)
% dim 2: HGF level
% dim 3: 1: muhat, 2: sahat, 3: mu, 4: sa

% Belief trajectories at 1st level

% x_state(r.irr) = [];
% y = r.y(:,1);
% y(r.irr) = [];

% Apply the unit-square sigmoid to the inferred states
prob = x_state.^ze./(x_state.^ze+(1-x_state).^ze);

% Initialize random number generator
if isnan(r.c_sim.seed)
    rng('shuffle');
else
    rng(r.c_sim.seed);
end

% Simulate responses
y_binary = binornd(1, prob);
yhat_binary = prob;


%% Run sim for continuous data modality (logRTs)

% mu1hat    = 1st level prediction
% sa1hat    = inverse precision of 1st level prediction
% mu2       = 2nd level outcome
% sa2       = inverse precisions of 2nd level outcome
% mu3       = 3rd level outcome


% Number of trials
n = size(infStates,1);

% Inputs
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
yhat = [yhat_binary yhat_reactionTime];


end

