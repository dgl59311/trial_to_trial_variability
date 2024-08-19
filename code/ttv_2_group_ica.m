%% Step 1: Group ICA
clc;
clear;
close all

%%
import_info

%% components for decomposition
ncomp_step_1    = 20;   % number of components to extract from individual data
ncomp_step_2    = 20;   % number of components to extract from stacked data
trial_length    = 1 * srate;  % trials are 1 second in total

%% find minimum number of trials available
available_trials = zeros(1, length(Participants));
parfor s = 1:size(Participants,1)
    available_trials(s) = count_trials(DataDir, Participants(s).name, []); 
    disp(s)
end

%% all subjects must have the same number of trials
% subjects 01_56_RV and s380 had less than 100 trials here they are removed
min_trials = 100;
Participants_m = Participants(find(available_trials >= min_trials));

%% tmp identify how much variance is explained with the 1st 20 components
tmpvariances = [];
parfor s = 1:size(Participants_m, 1)
    [eeg, ~, ~] = load_reshape(DataDir, Participants_m(s).name, min_trials, []);
    eeg = icatb_preproc_data(eeg, 'Remove Mean Per Timepoint');
    [coeff, score, explained, ~, ~, mu] = pca(eeg);
    exvar = (sum(explained(1:20))/sum(explained)) * 100;
    tmpvariances(s) = exvar;
end
% on average 20 components explain 87.9911 % of the variance

%% load times
load('times.mat')
index_1 = find(times >= 100, 1, 'first'); 
index_2 = find(times >= 300, 1,'first'); 
index_zero = find(times >= 0, 1, 'first');

%% loop through the participants and compute the PCA
firstPCA = [];
variances = [];
parfor s = 1:size(Participants_m, 1)
    [eeg, ~, ~] = load_reshape(DataDir, Participants_m(s).name, min_trials, []);
    eeg = icatb_preproc_data(eeg, 'Remove Mean Per Timepoint');
    % function for better performance obtained from eegift toolbox
    [whitesig_1, dewhiteM_1, Lambda, ~, ~] = icatb_calculate_pca(eeg', ...
        ncomp_step_1, 'whiten', true);   
    firstPCA(s).pca_.coeff = flip(dewhiteM_1'); 
    firstPCA(s).pca_.pcTime = flip(whitesig_1'); 
    latent_ = flip(diag(Lambda));
    variances(s, :) = cumsum(latent_)./sum(latent_);
    disp(s)
end

%%  build allData matrix for multi-level ICA/SOBI
allData = zeros(ncomp_step_1 * size(Participants_m, 1), min_trials * trial_length);
j = 1;
for i = 1:size(Participants_m,1)
   allData(j : ncomp_step_1 * i, :) = firstPCA(i).pca_.pcTime(1:ncomp_step_1,:);
   j = j + ncomp_step_1;
   disp(i) 
end

%% pca step 2 
allData = icatb_preproc_data(allData, 'Remove Mean Per Timepoint');
[whitesig_2, dewhiteM_2, ~, ~, ~] = icatb_calculate_pca(allData', ncomp_step_2,...
    'whiten', true);

%% perform mlGICA
[ica.weights, ica.sphere, ica.compvars, ica.bias, ica.signs, ica.lrates, ica.activations] = runica(whitesig_2', 'sphere', 'off');

%% obtain topographies from ICs
maps2 = inv(ica.weights)' * dewhiteM_2';
steps = 1:ncomp_step_2:size(dewhiteM_2, 1); 
% allocate data
topos = zeros(length(chanlocs), ncomp_step_2, length(Participants_m));
unmix = zeros(ncomp_step_2, length(chanlocs), length(Participants_m));

for i = 1:length(Participants_m) % loop through subjects
    topo = maps2(:,steps(i):steps(i)+ncomp_step_1-1) * firstPCA(i).pca_.coeff(1:ncomp_step_1,:);
    topos(:,:,i) = topo';
    unmix(:,:,i) = topo;
end

%% plot topographies
row = 5; 
col = 6;
comp = 1:ncomp_step_2;

figure()
for c = 1:size(comp,2)
    subplot(row,col,c)
    topoplot((mean(topos(:,comp(c),:),3)), chanlocs);
end

%% obtain activations
% compute group activations per condition
activity = reshape(ica.activations, ncomp_step_2, trial_length, min_trials);
activity = mean(activity,3);

%% plot IC-ERPs 
figure();
comp = 1:ncomp_step_2;
for c = 1:size(comp,2)
    subplot(5,4,c); hold on;
    plot(activity(c,:),'k');
    ylim([-4 4])
    axis tight;
end

%% compare IC with microstates
% load microstates data
prototype_maps = load(fullfile(ResultsFolder, '1_topographic_segmentations.mat'));
prototype_maps = prototype_maps.ALLEEG(9).microstate.prototypes(:, 1);

%% find the component that resemble the N1 microstate
top_IC = double(mean(topos,3));
cmaps = zeros(1, size(top_IC, 2));

for i = 1:size(top_IC, 2)    
    cmaps(i) = GlobalMapDiss(prototype_maps(:, 1), top_IC(:,i) * 1);
end
best_ic = find(cmaps == max(cmaps));

figure()
subplot(1, 2, 1)
topoplot(top_IC(:, best_ic), chanlocs);
subplot(1, 2, 2)
topoplot(prototype_maps(:, 1), chanlocs);

%% save data
times = get_timesvec(fullfile(DataDir, Participants(1).name));
save(fullfile(ResultsFolder, '2_mlGICA.mat'), 'ica', 'dewhiteM_2', 'topos', 'activity', 'unmix', 'times');

%% Step 2: Backfit to original data
% loop through the all the participants and backfit ICs
% try both cases, with group and single subject maps

%results_ = load(fullfile(ResultsFolder, '2_mlGICA.mat'));
%unmix = results_.unmix;

backfit_data = cell(length(Participants), 4);
fit_to_mean = 1;
var_map = var(mean(unmix, 3), [], 2);

for s = 1:size(Participants,1)
    % load data with all trials
    [eeg, n_trials, length_] = load_reshape(DataDir, Participants(s).name, [], []);
    
    if fit_to_mean == 1
        % fit the grand-average topographic map to the single trial data
        mean_um = mean(unmix, 3);
    else 
        % here it can also use the individual maps
        mean_um = unmix(:, :, s);
        for ii = 1:size(mean_um, 1)
            tmp_map = mean_um(ii, :);
            tmp_map = sqrt(var_map(ii) / var(tmp_map)) * tmp_map;
            % variance adjusted to mean variance of group map
            mean_um(ii, :) = tmp_map;
        end
    end
    
    backfit_ = mean_um * eeg;
    reshape_data = reshape(backfit_(best_ic, :), length_, n_trials);
    backfit_data{s, 1} = Participants(s).name;
    backfit_data{s, 2} = reshape_data;
    backfit_data{s, 3} = n_trials;
    backfit_data{s, 4} = length_;
    
end
%% save backfit data for alignment
save(fullfile(ResultsFolder, '2_backfit_maps'), 'backfit_data')

%%