function [logp, yhat, res] = obs_PerceptualBias_unitsq_sgm_tbt(r, infStates, ptrans)
% Calculates the log-probability of response y=1 under the unit-square sigmoid model
%
% --------------------------------------------------------------------------------------------------
% Copyright (C) 2012-2013 Christoph Mathys, TNU, UZH & ETHZ
%
% This file is part of the HGF toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.


% Predictions or posteriors?
pop = 1; % Default: predictions
if r.c_obs.predorpost == 2
    pop = 3; % Alternative: posteriors
end


% Transform zeta to its native space
ze = exp(ptrans(1));

% Initialize returned log-probabilities as NaNs so that NaN is
% returned for all irregualar trials
n = size(infStates,1);
logp = NaN(n,1);
yhat = NaN(n,1);
res  = NaN(n,1);

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
logp(reg) = y.*ze.*(logx -log1mx) +ze.*log1mx -log((1-x_state).^ze +x_state.^ze);
yhat(reg) = x_state;
res(reg) = (y-x_state)./sqrt(x_state.*(1-x_state));

return;
