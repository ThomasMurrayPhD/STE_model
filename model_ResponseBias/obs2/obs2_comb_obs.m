function [logp, yhat, res, logp_split] = obs2_comb_obs(r, infStates, ptrans)
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
ptrans_sgm = ptrans(1:2); % parameters for the sigmoid model
ptrans_logRT = ptrans(3:8); % parameters for the RT model


%% binary part of the response model

% compute log likelihood (binary responses)
[logp_binary, yhat_binary, res_binary] = ...
    obs2_unitsq_sgm_tbt(r, infStates, ptrans_sgm);

%% continuous part of the response model

% prepare inputs for logRT GLM assuming that the input is a matrix with 2
% columns (1: binary predictions, 2: log response times).
% rt = r;
% rt.y = r.y(:,2);

% Compute the log likelihood (logRTs)
[logp_reactionTime, yhat_reactionTime, res_reactionTime] = ...
    obs2_logrt_linear_binary(r, infStates, ptrans_logRT);


%% confidence part of the response model
% confidence01 = r.y(:,4);
% 
% [logp_confidence, yhat_confidence, res_confidence] = ...
%     tapas_softmax_binary(confidence01, infStates, ptrans);

%% get combined log likelihood of two response data modalities
logp = logp_binary + logp_reactionTime;
logp_split = [logp_binary logp_reactionTime];

%% return predictions and responses for each response data modality
yhat = [yhat_binary yhat_reactionTime];
res = [res_binary res_reactionTime];

end
