function [logp, yhat, res] = obs2_unitsq_sgm_tbt(r, infStates, ptrans)
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
zeta0 = exp(ptrans(1));
zeta1 = ptrans(2);

% Initialize returned log-probabilities as NaNs so that NaN is
% returned for all irregualar trials
n = size(infStates,1);
logp = NaN(n,1);
yhat = NaN(n,1);
res  = NaN(n,1);

% Weed irregular trials out from inferred states and responses
x_state = infStates(:,1,pop);

% add bias here
x_state(x_state == 0) = eps;
x_state(x_state == 1) = 1-eps;
x_state_temp = tapas_logit(x_state, 1);
x_state_temp = x_state_temp + (zeta1*r.u(:,2)); % add bias
x_state = tapas_sgm(x_state_temp, 1);


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
logp(reg) = y.*zeta0.*(logx -log1mx) +zeta0.*log1mx -log((1-x_state).^zeta0 +x_state.^zeta0);
yhat(reg) = x_state;
res(reg) = (y-x_state)./sqrt(x_state.*(1-x_state));

return;
