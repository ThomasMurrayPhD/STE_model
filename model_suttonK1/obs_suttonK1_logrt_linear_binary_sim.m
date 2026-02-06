function y = obs_suttonK1_logrt_linear_binary_sim(r, infStates, p)
% Simulates logRTs with Gaussian noise
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2016 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.

% Get parameters
be0  = p(1);
be1  = p(2);
be2  = p(3);
be3  = p(4);
be4  = p(5);
ze   = p(6);

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
y = logrt+sqrt(ze)*randn(n, 1);

end
