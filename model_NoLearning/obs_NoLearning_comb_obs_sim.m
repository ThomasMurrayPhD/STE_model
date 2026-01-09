function [y, yhat] = obs_NoLearning_comb_obs_sim(r, infStates, p)
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

% parameters for the psychometric model
alpha   = p(1); % PSE
beta    = p(2); % slope
gamma0  = p(3); % upper asymptote
lambda0 = p(4); % lower asymptote

% parameters for the gaussian RT model
mu      = p(5); % midpoint - JUST USE ALPHA
sigma0  = p(6); % width
gamma1  = p(7); % baseline - FIX IN CONFIG, DO NOT USE
lambda1 = p(8); % height
sigma1  = p(9); % noise


%% Run sim for binary predictions

% psychometric function
F = @(x,alpha,beta,gamma,lambda) gamma + (1 - gamma - lambda) ./ (1 + exp(-beta * (x - alpha)));

% inputs
u = r.u;

% loop through trials
yhat_binary = nan(size(u));
for i=1:numel(u)
    yhat_binary(i) = F(u(i), alpha, beta, gamma0, lambda0);
end

% Initialize random number generator
if isnan(r.c_sim.seed)
    rng('shuffle');
else
    rng(r.c_sim.seed);
end

% Simulate responses
y_binary = binornd(1, yhat_binary);



%% Run sim for continuous data modality (logRTs)

% gaussian function
G = @(x,mu,sigma0,gamma,lambda) gamma + lambda .* exp(-((x - mu).^2) ./ (2 * (sigma0.^2)));

% inputs
u = r.u;

% n
n = numel(u);

% Predict logRT
logrt=nan(size(u));
for i = 1:numel(u)
    % logrt(i) = G(u(i), mu, sigma0, gamma1, lambda1);
    logrt(i) = G(u(i), alpha, sigma0, gamma1, lambda1); % use alpha from psychometric function
end

% Initialize random number generator
if isnan(r.c_sim.seed)
    rng('shuffle');
else
    rng(r.c_sim.seed);
end

% Simulate (with noise)
y_RT = logrt+sqrt(sigma1)*randn(n, 1);
yhat_RT = logrt;




%% save values for both response data modalities
y = [y_binary y_RT];
yhat = [yhat_binary yhat_RT];

end

