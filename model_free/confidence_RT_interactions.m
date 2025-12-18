% script to run exploratory analyses with confidence/RT

% Get list of IDs
d = dir('..\STE_data\*.csv');
IDs = unique(cellfun(@(x) x(1:8), {d.name}, 'UniformOutput', false));

% Preallocate arrays
[safe_low, safe_high, threat_low, threat_high] = deal(nan(200,1));

% Loop through files
for i = 1:numel(d)
    tokens       = split(d(i).name, '_');
    ID           = tokens{1};
    condition    = extractBefore(tokens{3}, '.');

    data = readtable(['..\STE_data\', d(i).name]);


    if strcmp(condition, 'Safe')
        safe_low(strcmp(IDs, ID)) = mean(log(data.Response_RT(data.Confidence_idx==0)), 'omitnan');
        safe_high(strcmp(IDs, ID)) = mean(log(data.Response_RT(data.Confidence_idx==1)), 'omitnan');
    elseif strcmp(condition, 'Threat')
        threat_low(strcmp(IDs, ID)) = mean(log(data.Response_RT(data.Confidence_idx==0)), 'omitnan');
        threat_high(strcmp(IDs, ID)) = mean(log(data.Response_RT(data.Confidence_idx==1)), 'omitnan');
    end

end

%%
means = [mean(safe_low), mean(safe_high), mean(threat_low), mean(threat_high)];

figure; hold on;

scatter(ones(size(safe_low))    +((rand(size(safe_low))*0.5)-0.25),          safe_low);
scatter(ones(size(safe_high))   +((rand(size(safe_high))*0.5)-0.25) + 1,     safe_high);
scatter(ones(size(threat_low))  +((rand(size(threat_low))*0.5)-0.25) + 2,    threat_low);
scatter(ones(size(threat_high)) +((rand(size(threat_high))*0.5)-0.25) + 3,   threat_high);

%% anova
Y = [safe_low, safe_high, threat_low, threat_high];
T = array2table(Y, ...
    'VariableNames', {'SafeLow','SafeHigh','ThreatLow','ThreatHigh'});

Within = table( ...
    categorical({'Safe','Safe','Threat','Threat'})', ...
    categorical({'Low','High','Low','High'})', ...
    'VariableNames', {'Condition','Confidence'});



rm = fitrm(T, 'SafeLow-ThreatHigh ~ 1', 'WithinDesign', Within);

anovaTable = ranova(rm, 'WithinModel', 'Condition*Confidence');
disp(anovaTable)











