% script to plot which simulated parameters fail to recover


recov = importdata('model_suttonK1_recovery.mat');


failed = isnan(recov.LME);

param_names = fieldnames(recov);
param_names = param_names(1:end-3);

for i = 1:numel(param_names)
    figure('name', param_names{i}); hold on;
    S = recov.(param_names{i}).sim(~failed);
    F = recov.(param_names{i}).sim(failed);
    
    if strcmp(recov.(param_names{i}).space, 'log')
        S = log(S);
        F = log(F);
    end

    scatter(S, ones(size(S)) + .1, 'filled','blue');
    scatter(F, ones(size(F)), 'filled','red');
end



