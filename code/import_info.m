%% get data directories
DataDir = 'D:\trialvariability_sz_vbm_LOCAL\local_data';
ResultsFolder = 'D:\trialvariability_sz_vbm_LOCAL\local_results'; % use local folder
Participants = dir(fullfile(DataDir, '*.mat*'));
CurrDir = cd(); 
BehDir = 'C:\Users\gordillo\Dropbox\Code\trial_to_trial_variability_schizophrenia\preprocessing\behavior_eeg';

%% open eeglab
eeglabdir = 'D:\trialvariability_sz_vbm_LOCAL\software\eeglab2019_1';
cd(eeglabdir)
eeglab
cd(CurrDir)

%% read demographical files
demog_pat = readtable(fullfile(fileparts(CurrDir), ...
    'demographic_files', 'vbm_patients_2021_mf_v2.xlsx'));
demog_con = readtable(fullfile(fileparts(CurrDir), ...
    'demographic_files', 'vbm_controls_2021_mf_v2.xlsx'));

%% load time vector
load('times.mat')

%% set analysis
% read cap
ch_filename = fullfile(fileparts(CurrDir),...
    'preprocessing', '64-4_Biosemi.xyz');    % cap coordinates filename
chanlocs = readlocs(ch_filename, 'filetype', 'xyz');    % reads the channels location with EEGLAB
chanlocs(65:68) = [];   % remove ocular channels

%% info of dataset
srate = 512;
nbchan = 64;