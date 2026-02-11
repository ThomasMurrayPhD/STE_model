function [y, logrt] = obs_PredictionBias_logrt_linear_binary_sim(r, infStates, p)
% [y, logrt] = m1_logrt_linear_binary_sim(r, infStates, p)
%
% Simulates logRTs with Gaussian noise.
% (Designed to be compatible with the HGF Toolbox as part of TAPAS).
%
% INPUT
%   r             struct      Struct obtained from tapas_simModel.m
%   infStates     tensor      Tensor containing inferred states from the
%                             perceptual model    
%   p             vector      1xP vector with free param values (nat space)
%
%   OPTIONAL:
%
% OUTPUT    
%   y             vector       Nx1 vector with simulated logRTs
%   logrt         vector       Nx1 vector containing noise-free predictions
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





% mu1hat    = 1st level prediction
% sa1hat    = inverse precision of 1st level prediction
% mu2       = 2nd level outcome
% sa2       = inverse precisions of 2nd level outcome
% mu3       = 3rd level outcome



% Get parameters
be0  = p(1);
be1  = p(2);
be2  = p(3);
be3  = p(4);
be4  = p(5);
sa   = p(6);

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


% mu1 = infStates(:,1,3);
% pu = 0.5-abs(mu1-0.5);


% Calculate predicted log-reaction time
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% logrt = be0 +be1.*surp +be2.*bernv +be3.*inferv +be4.*pv;
logrt = be0 +be1.*surp +be2.*bernv +be3.*inferv +be4.*stim_noise;

% logrt = be0 +be1.*surp +be2.*bernv +be3.*inferv +be4.*pu; %%% THIS SEEMS
% TO WORK... MAX RT SHIFTS WITH PSE


% Initialize random number generator
if isnan(r.c_sim.seed)
    rng('shuffle');
else
    rng(r.c_sim.seed);
end

% Simulate
y = logrt+sqrt(sa)*randn(n, 1);

end
