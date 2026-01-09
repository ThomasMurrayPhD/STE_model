function [logp, yhat, res, logp_split] = obs_NoLearning_comb_obs(r, infStates, ptrans)
% [logp, yhat, res, logp_split] = obs3_comb_obs(r, infStates, ptrans)
%
% Calculates the combined log-probability of binary and continuous
% behavioural responses (NO LEARNING)
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


%% Separate parameters


% Transform parameters into native space
% parameters for the psychometric model
alpha   = tapas_sgm(ptrans(1), 1); % PSE
beta    = exp(ptrans(2)); % slope
gamma0  = tapas_sgm(ptrans(3), 1); % upper asymptote
lambda0 = tapas_sgm(ptrans(4), 1); % lower asymptote

% parameters for the gaussianRT model
mu      = tapas_sgm(ptrans(5), 1); % midpoint - JUST USE ALPHA
sigma0  = exp(ptrans(6)); % width
gamma1  = exp(ptrans(7)); % baseline - FIX IN CONFIG, DO NOT USE
lambda1 = exp(ptrans(8)); % height
sigma1  = exp(ptrans(9)); % noise


%% binary part of the response model

% psychometric function
F = @(x,alpha,beta,gamma,lambda) gamma + (1 - gamma - lambda) ./ (1 + exp(-beta * (x - alpha)));

% inputs
u = r.u;

% responses
y = r.y(:,1);

% probabilities
yhat_binary=nan(size(u));
logp_binary=nan(size(u));

for i = 1:numel(u)
    % p(response==1)
    p_1 = F(u(i),alpha,beta,gamma0,lambda0);

    % predicted outcome
    yhat_binary(i)=p_1;

    % probability of observed outcome
    logp_binary(i) = r.y(i) * log(p_1) + (1 - r.y(i)) * log(1 - p_1);

end

% Calculate residuals (standardized difference between observed and predicted)
res_binary = (y - yhat_binary) ./ sqrt(yhat_binary .* (1 - yhat_binary));




%% continuous part of the response model

% gaussian function
G = @(x,mu,sigma,gamma,lambda)gamma+(lambda*exp(-((x-mu)^2)/(2*(sigma^2))));

% inputs
u = r.u;

% responses
y = r.y(:,2);

% n
n = numel(u);

% probabilities
[yhat_RT, logp_RT, res_RT] = deal(nan(size(u)));

% Predict logRT
logrt=nan(size(u));
for i = 1:n
    % logrt(i) = G(u(i), mu, sigma0, gamma1, lambda1);
    logrt(i) = G(u(i), alpha, sigma0, gamma1, lambda1); % use alpha from psychometric function
end

% Calculate log-probabilities (taken straight from original tapas model)
% 
% Note: 8*atan(1) == 2*pi (this is used to guard against
% errors resulting from having used pi as a variable).
reg = ~ismember(1:n,r.irr);
logp_RT(reg) = -1/2.*log(8*atan(1).*sigma1) -(y-logrt).^2./(2.*sigma1);
yhat_RT(reg) = logrt;
res_RT(reg) = y-logrt;


%% confidence part of the response model




%% get combined log likelihood of two response data modalities
logp = logp_binary + logp_RT;
logp_split = [logp_binary logp_RT];

%% return predictions and responses for each response data modality
yhat = [yhat_binary yhat_RT];
res = [res_binary res_RT];

end
