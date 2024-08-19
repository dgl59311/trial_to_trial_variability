%% step 3: align n1 peak
clc;
clear;
close all;

%%
import_info

%% alignment (these values are taken from the tutorial of the toolbox) 
options.verbose = false; % show figures
options.disp_log = false; % print log messages
options.show_pca = false;
options.use_maximum = false;
options.sigma = [0.01:0.01:0.2];
options.alpha = [0.001,0.01,0.1];

%% time window where to find the peak
time_0 = 1;
time_1 = length(times);
index_1 = find(times >= 100, 1, 'first'); 
index_2 = find(times >= 300, 1,'first'); 

%% get id of subjects for results
id_subjects = {};
for i = 1:length(Participants)
    id_ = Participants(i).name;
	ix_d = strfind(id_, '_pat_');
    if isempty(ix_d)
        ix_d = strfind(id_, '_con_');
    end
	id_ = strrep(id_, id_(ix_d:end), '');
    id_subjects{i} = id_;
end

clear id
id_subjects = id_subjects';

%% loop starts
% use group ic maps
load(fullfile(ResultsFolder, '2_backfit_maps.mat'))    
ResultsMatrix = zeros(length(Participants), 13);

%%
tic
parfor j = 1:length(Participants)

    id_ = Participants(j).name;
    data_ = backfit_data{j, 2};
    data_ = data_ * 1;
    eegset = set_eegalignment(data_, chanlocs,...
        times, size(data_, 2), srate);
    [EEG, com, order, lags, event_type, E_lags] = pop_extractlag(eegset, 0, 1,...
        [times(index_1), times(index_2)], options);

    % ms to sample
    [lag_sample] = ttv_ms_to_sample(lags, times);

    % finding of minimum values + making cell array of all 
    min_values = zeros(1, length(lags));
    for k = 1 : length(lags)
        min_values(k) = data_(lag_sample(k), k);
    end 

    ntrials = length(min_values);
    % variables to retain 
    unaligned = mean(data_, 2);
    non_aligned = min(unaligned(index_1:index_2));
    aligned = mean(min_values);
    cv = (std(lags)/mean(lags)) * 100;
    lag = mean(lags);

    % robust measures of variability
    % quartile CV
    Q1 = prctile(lags, 25);
    Q3 = prctile(lags, 75);
    qcv = (Q3 - Q1) / (Q3 + Q1);
    % robust cv
    rcv = ((Q3 - Q1) / median(lags));

    % indicate condition

    if contains(id_, '_l_soa')
        condition = 1;
    elseif contains(id_, '_s_soa')
        condition = 2;
    elseif contains(id_, '_mask')
        condition = 3;
    elseif contains(id_, '_vernier')
        condition = 4;
    end      

    % write table 
    if contains(id_, 'pat')
        ix_d = strfind(id_, '_pat_');
        idmf_ = strrep(id_, id_(ix_d:end), '');
        tab = demog_pat(contains(demog_pat.matName, idmf_), :);  
        group = 1;
    else 
        ix_d = strfind(id_, '_con_');
        idmf_ = strrep(id_, id_(ix_d:end), '');
        tab = demog_con(contains(demog_con.matName, idmf_), :);
        group = 2;
    end

    gender = tab.Gender;
    age = tab.Age;
    education = tab.Education;

    store_ = [NaN; group; gender; age; education; condition; aligned; ...
        non_aligned; cv; lag; ntrials; rcv; qcv]';
    ResultsMatrix(j, :) = store_;

end
toc

%% write csv table
main_results = array2table(ResultsMatrix);
main_results.ResultsMatrix1 = id_subjects;
main_results.Properties.VariableNames = {'ID', ...
    'Group', 'Gender', 'Age', 'Education', 'Condition', 'MinAl', 'MinNonAl',...
    'CV', 'Lags', 'NTrials', 'RobustCV', 'QCV'};
writetable(main_results,...
    fullfile(ResultsFolder, ['3_results_ica_alignment.csv']));

%%