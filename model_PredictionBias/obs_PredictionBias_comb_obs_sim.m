function [y, yhat] = obs1_comb_obs_sim(r, infStates, p)
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

p_sgm = p(1); % parameter for the sigmoid model
p_logRT = p(2:7); % parameters for the RT model


%% Run sim for binary predictions
[pred, yhat_pred] = obs1_unitsq_sgm_tbt_sim(r, infStates, p_sgm);

%% Run sim for continuous data modality (logRTs)
[logReactionTime, yhat_rt] = obs1_logrt_linear_binary_sim(r, infStates, p_logRT);


%% save values for both response data modalities
y = [pred logReactionTime];
yhat = [yhat_pred yhat_rt];

end

