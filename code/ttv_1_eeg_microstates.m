%% topographic segmentation analysis
clc;
clear;
close all;

%%
import_info

%% 
% Separate the participants by group 
groups = {'_pat_', '_con_'};
[split_data] = find_in_cell(Participants, groups);
all_patients = split_data{1};
all_controls = split_data{2};

%% obtain grand means by condition
% provide a directory where behavior data are stored if you want to
% consider only hit trials; here we use all trials
% Controls
[EEGFiles{1}] = get_grand_means(all_controls, {'controls'}, {'_vernier'}, []);
[EEGFiles{2}] = get_grand_means(all_controls, {'controls'}, {'_l_soa'}, []);
[EEGFiles{3}] = get_grand_means(all_controls, {'controls'}, {'_s_soa'}, []);
[EEGFiles{4}] = get_grand_means(all_controls, {'controls'}, {'_mask'}, []);

% Patients
[EEGFiles{5}] = get_grand_means(all_patients, {'patients'}, {'_vernier'}, []);
[EEGFiles{6}] = get_grand_means(all_patients, {'patients'}, {'_l_soa'}, []);
[EEGFiles{7}] = get_grand_means(all_patients, {'patients'}, {'_s_soa'}, []);
[EEGFiles{8}] = get_grand_means(all_patients, {'patients'}, {'_mask'}, []);

%% concatenate files for segmentation
ALLEEG = [];
for i=1:length(EEGFiles)
    g_avg = delete_bl(EEGFiles{i});
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, g_avg, 0 );
    eeglab redraw % updates EEGLAB datasets
end
 
%% select data for microstates analysis
[EEG, ALLEEG] = pop_micro_selectdata( EEG, ALLEEG, 'datatype', 'ERPavg',... 
'avgref', 1, ...
'normalise', 1, ...
'dataset_idx', 1:length(EEGFiles) );
% store data in a new EEG structure
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw % updates EEGLAB datasets
 
%% microstate segmentation
% Select the "GFPpeak" dataset and make it the active set
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG,length(EEGFiles) ,'retrieve',length(EEGFiles)+1,'study',0);
eeglab redraw
 
% Perform the microstate segmentation
EEG = pop_micro_segment( EEG, 'algorithm', 'kmeans', ...
'Nmicrostates', 2:8);
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
 
%%
figure; 
MicroPlotTopo( EEG, 'plot_range', [] );
%%
EEG = pop_micro_selectNmicro( EEG );
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% backfit data
for i = 1:length(EEGFiles)
    
    fprintf('Importing prototypes and backfitting for dataset %i\n',i)
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',i,'study',0);
    EEG = pop_micro_import_proto( EEG, ALLEEG, length(EEGFiles)+1);
     
    % Back-fit microstates on EEG
    EEG = pop_micro_fit( EEG, 'polarity', 1 );
    % Temporally smooth microstates labels
    EEG = pop_micro_smooth( EEG, 'label_type', 'backfit',... 
    'smooth_type', 'reject segments', ...
    'minTime', 15, ...
    'polarity', 1 );
    % Calculate microstate statistics
    EEG = pop_micro_stats( EEG, 'label_type', 'backfit',... 
    'polarity', 1 );
	[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
end

%%
save(fullfile(ResultsFolder, '1_topographic_segmentations.mat'), 'ALLEEG')

%% plotting GFP of active microstates for the first 1500 ms for each condition.
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',1,'study',0);
figure;
MicroPlotSegments( EEG, 'label_type', 'backfit', ...
'plotsegnos', 'first', 'plot_time', [0 490], 'plottopos', 1);
eeglab redraw
 
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',2,'study',0);
figure;
MicroPlotSegments( EEG, 'label_type', 'backfit', ...
'plotsegnos', 'first', 'plot_time', [0 490], 'plottopos', 1);
eeglab redraw
 
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',3,'study',0);
figure;
MicroPlotSegments( EEG, 'label_type', 'backfit', ...
'plotsegnos', 'first', 'plot_time', [0 490], 'plottopos', 1);
eeglab redraw
 
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',4,'study',0);
figure;
MicroPlotSegments( EEG, 'label_type', 'backfit', ...
'plotsegnos', 'first', 'plot_time', [0 490], 'plottopos', 1);
eeglab redraw
 
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',5,'study',0);
figure;
MicroPlotSegments( EEG, 'label_type', 'backfit', ...
'plotsegnos', 'first', 'plot_time', [0 490], 'plottopos', 1);
eeglab redraw
 
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',6,'study',0);
figure;
MicroPlotSegments( EEG, 'label_type', 'backfit', ...
'plotsegnos', 'first', 'plot_time', [0 490], 'plottopos', 1);
eeglab redraw
 
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',7,'study',0);
figure;
MicroPlotSegments( EEG, 'label_type', 'backfit', ...
'plotsegnos', 'first', 'plot_time', [0 490], 'plottopos', 1);
eeglab redraw

[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',8,'study',0);
figure;
MicroPlotSegments( EEG, 'label_type', 'backfit', ...
'plotsegnos', 'first', 'plot_time', [0 490], 'plottopos', 1);
eeglab redraw

%%



