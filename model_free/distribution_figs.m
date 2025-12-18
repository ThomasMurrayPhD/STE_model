% Script to visualise psychometric, RT and confidence distributions

% Get list of IDs
d = dir('..\STE_data\*.csv');
IDs = unique(cellfun(@(x) x(1:8), {d.name}, 'UniformOutput', false));

% Preallocate arrays
data.safe.psad      = deal(nan(numel(IDs), 6));
data.safe.rt        = deal(nan(numel(IDs), 6));
data.safe.conf      = deal(nan(numel(IDs), 6));
data.threat.psad    = deal(nan(numel(IDs), 6));
data.threat.rt      = deal(nan(numel(IDs), 6));
data.threat.conf    = deal(nan(numel(IDs), 6));


% Loop through files
for i = 1:numel(d)
    tokens       = split(d(i).name, '_');
    ID           = tokens{1};
    condition    = extractBefore(tokens{3}, '.');

    sub_data = readtable(['..\STE_data\', d(i).name]);

    % get p(sad)
    psad = arrayfun(@(x) mean(sub_data.Response_idx(sub_data.Outcome_p_sad==x), 'omitnan'), 0:20:100);

    % get RT
    logrt = log(arrayfun(@(x) mean(sub_data.Response_RT(sub_data.Outcome_p_sad==x), 'omitnan'), 0:20:100));

    % Get conf
    conf = arrayfun(@(x) mean(sub_data.Confidence_idx(sub_data.Outcome_p_sad==x), 'omitnan'), 0:20:100);

    % Add to matrices
    if strcmp(condition, 'Safe')
        data.safe.psad(strcmp(IDs, ID), :)      = psad;
        data.safe.rt(strcmp(IDs, ID), :)        = logrt;
        data.safe.conf(strcmp(IDs, ID), :)      = conf;
    elseif strcmp(condition, 'Threat')
        data.threat.psad(strcmp(IDs, ID), :)    = psad;
        data.threat.rt(strcmp(IDs, ID), :)      = logrt;
        data.threat.conf(strcmp(IDs, ID), :)    = conf;
    end

end


%% Visualise

conds = {'safe', 'threat'};
vars = {'psad', 'rt', 'conf'};

for iC = 1:2
    figure('name', conds{iC});
    
    for iV = 1:3
        subplot(3,1,iV);
        hold on;
        
        % get relevant data
        d = data.(conds{iC}).(vars{iV});

        % calculate mean
        M = mean(d, 1, 'omitnan');
        
        % plot data with jitter
        for i=1:size(d,1)
            scatter(...
                arrayfun(@(x) x+((rand-0.5)/3), 1:6), ...
                d(i,:), 5, [.7,.7,.7]);
        end
        
        % plot mean
        plot(1:6, M, '-o',...
            'color', [0,0.4470,0.7410],...
            'linewidth', 1.5, 'MarkerSize',8,...
            'MarkerEdgeColor',[0,0.4470,0.7410],...
            'MarkerFaceColor',[1,1,1])
        
        % set properties
        set(gca, 'xlim', [1,6], 'xtick', 1:6, 'xticklabels', {'0', '20', '40', '60', '80', '100'});

    end
end


%% Plot safe/threat difference

figure('name', 'safe-threat difference');


for iV = 1:3
    subplot(3,1,iV);
    hold on;
    
    % get safe and threat data
    s = data.safe.(vars{iV});
    t = data.threat.(vars{iV});

    d = s-t;
    
    % calculate mean
    M = mean(d, 1, 'omitnan');
    
    % plot data with jitter
    for i=1:size(d,1)
        scatter(...
            arrayfun(@(x) x+((rand-0.5)/3), 1:6), ...
            d(i,:), 5, [.7,.7,.7]);
    end
    
    % plot mean
    plot(1:6, M, '-o',...
        'color', [0,0.4470,0.7410],...
        'linewidth', 1.5, 'MarkerSize',8,...
        'MarkerEdgeColor',[0,0.4470,0.7410],...
        'MarkerFaceColor',[1,1,1])
    
    % set properties
    set(gca, 'xlim', [1,6], 'xtick', 1:6, 'xticklabels', {'0', '20', '40', '60', '80', '100'});

end












