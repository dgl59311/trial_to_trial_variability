%%APP for Trial-to-trial variability study V2. (increase epoch sizes)
% based on the adapted version by Simona:APP_BP_Oph_SG.m
% eeglab to be installed and set to path
clc; 
clear; 
close all;
CurrDir = pwd;    % gets current directory
% Paths to the data
ParticipantsDir = uigetdir([],'Path to the raw data of all group of participants of the experiment');
% Path to save data
SaveDir   = uigetdir([],'Path to store the results');

%% start eeglab

%%
% Info about your data
n_ch = 64;                              % number of EEG channels
n_runs = 8;                             % number of runs
n_cond = 4;                             % number of conditions
%events = num2cell(100:199);            % Stimulus events of interest
epoch = [-0.5 0.5];                         % boundaries of the epoch related to stimulus onset
prestim = [epoch(1) 0];                 % pre-stimulus to be removed
eog_th = 100;                           % eog threshold around 75-100 uV

% Info about the analysis
doICA = 1;                              % if want to do ICA decomposition, set doICA = 1
doDownSamp = 1;                         % if want to do Downsampling, set doDownSamp = 1
d_Fs = 512;                             % desired frequency in Hz
doFiltering = 1;                        % if want to do Filtering, set doFiltering = 1
LowCutOffFreq = 1;                      % lower bound cut-off frequency
HighCutOffFreq = 40;                    % high bound cut-off frequency
filterOrder = [];                       % filterOrder - leave blank and determined by eeglab
doRemPrestBaseline = 1;                 % if wants to remove prestimulus baseline, set prest = 1

% EOG cap
CapDir = CurrDir;                       			% directory with the caps in xyz coordinates % CapDir = 'D:\3D caps\EEGLAB';
cd(CapDir);                             			% changes to the directory with the cap
ch_filename = '64-4_Biosemi.xyz';       			% cap coordinates filename
chanlocs = readlocs(ch_filename,'filetype','xyz'); 	% reads the channels location with EEGLAB
cd(CurrDir);
eye_conv = 2;                           			% if Lab EOG convention eye_conv = 1; if Georgia = 2

%% 
Group={'patients','controls'}';

for iGroup = 1:size(Group, 1)
    
    % Group directory with the data
    GroupDir = strcat(fullfile(ParticipantsDir,char(Group{iGroup})));
    
    % Creates the directory to save data for each group
    mkdir(SaveDir,Group{iGroup});
    
    % Subject folders
    SubjFolders = dir(GroupDir);
    SubjFolders(1:2) = [];
    isub = [SubjFolders(:).isdir];
    % The pool of subjects
    Subject_pool = {SubjFolders(isub).name}';
    
    
    for iSubject = 1:size(Subject_pool,1)
        
        % Subject directory
        SubjDir = strcat(GroupDir,'\',char(Subject_pool{iSubject}));
        cd(SubjDir);
        % BDF & SING folders
        BDFFolder = dir('*orig*');
        SINGFolder = dir('*SING*');
        
       try
            % change to BDF folder and get the name of "run" EEG BDF files
            cd(strcat(SubjDir,'\',BDFFolder.name));
            bdf = dir('*run*');
            bdfName = {bdf(:).name}';
            
            % assert the runs are in order because there are some inconsistency
            % in names in the Georgia files
            num_run = cellfun(@(x) x(end-4), bdfName, 'UniformOutput',false);
            [~,ind_sort] = sort(str2num(cell2mat(num_run)));
            bdfName = {bdfName{ind_sort}}';
            
            % changes to SING folder and get the name of "run" DV files
            cd(strcat(SubjDir,'\',SINGFolder.name));
            dv = dir('*dv*');
            dvName = {dv(:).name}';
            % assert the runs are in order
            num_run = cellfun(@(x) x(end-3), dvName, 'UniformOutput',false);
            [~,ind_sort] = sort(str2num(cell2mat(num_run)));
            dvName = {dvName{ind_sort}}';
            
            % check if the dv, bdf and number of runs are the same
            % sometimes some dataset does not have a bdf or a dv file
            if size(bdfName,1) ~= n_runs || size(dvName,1) ~= n_runs
                disp('The number of dv and bdf files are not correspondent to the number of runs');
                disp(strcat('Problem with subject: ',char(Subject_pool{iSubject})));
                cd(SaveDir);
                save(strcat('error_',char(Subject_pool{iSubject})),'bdfName','dvName');
                % Intersect the dv and bdf files and get the ones that are
                % corresponding and continue the analysis only those are in the
                % DV and BDF
                dv_run = cellfun(@(x) x(end-3), dvName, 'UniformOutput',false);
                bdf_run = cellfun(@(x) x(end-4), bdfName, 'UniformOutput',false);
                [~,dv_run_in,bdf_run_in] = intersect(dv_run,bdf_run);
                dvName = {dvName{dv_run_in}}';
                bdfName = {bdfName{bdf_run_in}}';
            end
            
            % find where the bdf file has a discontinuity in the triggers
            error_trials = 0;
            % Bad channels
            bad_channels = 0;
            % Invalid trials
            invalid_trials = 0;
            % Bad trials
            bad_trials = 0;
            % Blink trials
            blink_trials = 0;
            % Bad channels in each trial
            trial_bad_ch = 0;
            
            
            for iRun = 1:n_runs;
                
                % name of the dv and bdf files to be read
                dvFile = dvName{iRun};
                bdfFile = bdfName{iRun};
                
                % changes to the SING folder and reads the dv file
                cd(strcat(SubjDir,'\',SINGFolder.name));
                dv = lpsy.readDvFile(dvFile);               % reads invalid trials as well
                data_dv = dv.pool0;
                % the lab convention is from 100 to 199; however, sometimes
                % some subjects have too many invalid trials and we can exceed
                % 199. Since the trial counts start back from 100, this plus199
                % function takes care of it
                cd(CurrDir);
                data_dv(:).trigger_1 = plus199(data_dv(:).trigger_1);
                data_dv_events = data_dv(:).trigger_1; % all trial beginnings events in dv file
                
                % changes to the BDF folder and reads the bdf file using eeglab
                cd(strcat(SubjDir,'\',BDFFolder.name));
                % Read in the BDF file using eeglab
                dataset = pop_biosig(bdfFile, 'channels', 1: n_ch+4); % n eeg channels and 4 EOG
                
                % Allocates the channel locations to the EEG structure
                dataset.chanlocs = chanlocs;
                
                % Check if the recording ended too early or too late
                % Std is in the range 100~199
                dataset_events = [dataset.event(:).type]';
                %         [dataset_events, error_trial] = qualityCheck(dataset_events);  % the error trial will later be removed from eeg
                %         error_trials = error_trials + length(error_trial);
                Ind = find(dataset_events>=100 & dataset_events<=250); % new, experimental
                dataset_events = dataset_events(Ind);
                cd(CurrDir);
                dataset_events = plus199(dataset_events); % new, this is experimental

                j=1;
                for i=[Ind]' % new
                    dataset.event(i).type = dataset_events(j);
                    dataset.urevent(i).type = dataset_events(j);
                    j=j+1;
                end
                
                % check if EEG and DV correspond
                % intersect both and remove the trials that are not on both
                [dataset_events,data_dv_in,dataset_events_in] = intersect(data_dv_events,dataset_events);
                for f=fieldnames(data_dv)'
                    data_dv.(f{1})=data_dv.(f{1})(data_dv_in);
                end
				% Here the pre-processing starts
				% ------------------------------------------------------------------------------------------------------------------------------------------
                % Bandpass filter with 0 phase lag
                if doFiltering == 1
				
                    dataset = pop_eegfiltnew(dataset, LowCutOffFreq, HighCutOffFreq, filterOrder, 0, [], 0);
                    
                    % CleanLine to remove the 50 Hz and the harmonics
                    %dataset = pop_cleanline(dataset, 'Bandwidth',2,'ChanCompIndices',[1:dataset.nbchan],...
                    %    'SignalType','Channels','ComputeSpectralPower',true,               ...
                    %    'LineFrequencies',[50 100 150 200 250] ,'NormalizeSpectrum',false, ...
                    %    'LineAlpha',0.01,'PaddingFactor',2,'PlotFigures',false,            ...
                    %    'ScanForLines',true,'SmoothingFactor',100,'VerboseOutput',1,       ...
                    %    'SlidingWinLength',4,'SlidingWinStep',4); close all;
                end
                
                % Downsample the data if needed
                if doDownSamp == 1 && dataset.srate ~= d_Fs
                    dataset = pop_resample(dataset,d_Fs);
                end
                
                % Calculate the EOG channels, increase the SNR
                if eye_conv == 1
                    % using the Lab convention
                    veog = dataset.data(n_ch+1,:) - dataset.data(n_ch+2,:); % VEOG channel (top - bottom)
                    heog = dataset.data(n_ch+3,:) - dataset.data(n_ch+4,:); % HEOG channel (left-right - looking at the subj)
                elseif eye_conv == 2
                    % using the Georgia convention
                    heog = dataset.data(n_ch+2,:) - dataset.data(n_ch+1,:); % HEOG channel (left-right - looking at the subj)
                    veog = dataset.data(n_ch+3,:) - dataset.data(n_ch+4,:); % VEOG channel (top - bottom)
                end
                
                dataset.data(n_ch+1,:) = veog; % overwrite channels
                dataset.data(n_ch+2,:) = heog;
                clear veog heog % free some memory
                
                dataset = pop_select(dataset, 'channel', 1:n_ch+2); % remove channels 67 68
                %rename the EOG channels
                dataset.chanlocs(n_ch+1).labels = 'VEOG';
                dataset.chanlocs(n_ch+2).labels = 'HEOG';
                % extract EOG channels
                dataset_eog = pop_select(dataset, 'channel', n_ch+1:n_ch+2);
                % remove EOG channels from EEG data
                dataset = pop_select(dataset, 'channel', 1:n_ch);
                
                %------------------------------------------------------------------------------------
                
                % Re-reference to the biweight mean of the channels
                % Because according to Biosemi the CMS electrode does not provide
                % the full 80 dB CMRR
                [avg_ch,~] = myBiweight(dataset.data');
                dataset.data = dataset.data - repmat(avg_ch,n_ch,1);
                clear avg_ch
                
                %------------------------------------------------------------------------------------
                
                % Find bad channels
                % standard deviation
                [~,ch_std] = myBiweight(dataset.data);
                bad_ch_ind_1 = myFindOutliers(ch_std);
                % mean correlation of top 4 correlation coefficients
                ch_r_raw = corr(dataset.data');
                ch_r_raw = sort(ch_r_raw);
                ch_r = mean(ch_r_raw(end-4:end-1,:)); % top 4 correlation coefficients excluding self-correlation
                bad_ch_ind_2 = myFindOutliers(ch_r);
                % remove and interpolate bad channels
                bad_ch = unique([bad_ch_ind_1,bad_ch_ind_2]); % indeces of bad channels
                if ~isempty(bad_ch)
                    dataset = eeg_interp(dataset, bad_ch, 'spherical');
                end
                bad_channels = bad_channels + length(bad_ch); % the number of bad channels per subject
                %------------------------------------------------------------------------------------
                
                % Epoch data
                [dataset, indices_d] = pop_epoch( dataset, num2cell(dataset_events), epoch); % convention: events 100 to 199 mark stimulus onset changed, no upper limit
                dataset = eeg_checkset( dataset );
                [dataset_eog, indices] = pop_epoch( dataset_eog, num2cell(dataset_events), epoch);
                dataset_eog = eeg_checkset( dataset_eog );  
                % NEW: if when epoching, data are deleted, also delete from
                % dv file
                if length(indices_d) ~= length(dataset_events)
                    missingvals = setdiff(1:length(dataset_events), indices_d);
                    for f=fieldnames(data_dv)'
                        data_dv.(f{1})(missingvals)=[];
                    end
                end
                
                %----------------------------------------------------------------------------------
                % NECESSARY THIS SECTION???
                
                % remove invalid epochs
                inval_epoch1 = find(data_dv.valid == 0); % find the invalid epochs
                inval_epoch2 = find(data_dv.react_ti < 300 | data_dv.react_ti > 3000); % get react_ti too long or too short
                inval_epoch = unique([inval_epoch1' inval_epoch2']);
                dataset = pop_select(dataset,'notrial',inval_epoch);
                dataset_eog = pop_select(dataset_eog,'notrial',inval_epoch);
                clear inval_epoch1 inval_epoch2 % free some memory
                invalid_trials = invalid_trials + length(inval_epoch);
                % remove the invalid trials from the dv data
                for f=fieldnames(data_dv)'
                    data_dv.(f{1})(inval_epoch)=[];
                end
                
                % remove trials where the subject did not see the stimulus
                % the eog during stimulus presentation is larger than 75uV
                % find point time equal to -100 and 100 ms since probably the most
                % important value
                [~,t_100] = min(abs(dataset_eog.times+100)); % - 100 ms
                [~,t100] = min(abs(dataset_eog.times-100)); % 100 ms
                eyeblink_trial = [];
                for itrial = 1:dataset_eog.trials
                    v_blink = find(dataset_eog.data(1,t_100:t100,itrial) > eog_th); % vertical blink
                    h_blink = find(dataset_eog.data(2,t_100:t100,itrial) > eog_th); % saccade
                    if ~isempty(v_blink) | ~isempty(h_blink)
                        eyeblink_trial = [eyeblink_trial itrial];
                    end
                end
                % remove the epochs of blink trials
                dataset = pop_select(dataset, 'notrial', eyeblink_trial);
                dataset_eog = pop_select(dataset_eog, 'notrial', eyeblink_trial);
                % remove the dv data of blink trials
                for f=fieldnames(data_dv)', data_dv.(f{1})(eyeblink_trial)=[]; end
                
                blink_trials = blink_trials + length(eyeblink_trial);
                clear eyeblink_trial % free some memory
                
                %----------------------------------------------------------------------------------
                
                % Remove trials that might contain artifacts that are too big for ICA decomposition
                % Remember the "garbage in, garbage out"
                % Hopefully not many :)
                
                % Max amplitude difference
                amp_diffs = zeros(size(dataset.data,1),size(dataset.data,3));
                for iChan = 1:size(dataset.data,1)
                    for itrial = 1:size(dataset.data,3)
                        amp_diffs(iChan,itrial) = max(dataset.data(iChan,:,itrial)) - min(dataset.data(iChan,:,itrial));
                    end
                end
                [epoch_amp_d,~] = myBiweight(amp_diffs');
                % Epoch variance or the mean GFP
                epoch_GFP = mean(squeeze(std(dataset.data,0,2)));
                % Epoch's mean deviation from channel means.
                [means,~] = myBiweight(dataset.data(:,:)); % channel mean for all epochs
                epoch_m_dev = zeros(1,size(dataset.data,3));
                for itrial = 1:size(dataset.data,3)
                    epoch_m_dev(itrial) = mean(abs(squeeze(mean(dataset.data(:,:,itrial),2))' - means));
                end
                
                % Find the bad trials
                Rej_ep_amp_d = myFindOutliers(epoch_amp_d);
                Rej_ep_GFP = myFindOutliers(epoch_GFP);
                Rej_ep_mdev = myFindOutliers(epoch_m_dev);
                Rej_epoch = unique([Rej_ep_amp_d Rej_ep_GFP Rej_ep_mdev]); % indeces of epochs rejected
                
                % Remove the bad trials
                dataset = pop_select(dataset,'notrial',Rej_epoch);
                dataset_eog = pop_select(dataset_eog,'notrial',Rej_epoch);
                % Count the bad trials
                bad_trials = bad_trials + length(Rej_epoch);
                
                
                % Add epochs to corresponding conditions
                for f=fieldnames(data_dv)'
                    data_dv.(f{1})(Rej_epoch)=[];
                end
                
                
                if doICA == 1
                    %------------------------------------------------------------------------------------
                    % Perform ICA - SOBI
                    %------------------------------------------------------------------------------------
                    EEG = pop_runica(dataset, 'icatype', 'sobi', 'dataset',1, 'options',{});
                    if isempty(EEG.icaact)
                        disp('EEG.icaact not present. Recomputed from data.');
                        if length(size(EEG.data))==3
                            EEG.icaact = reshape(EEG.icaweights*EEG.icasphere*reshape(EEG.data,[size(EEG.data,1)...
                                size(EEG.data,2)*size(EEG.data,3)]),[size(EEG.data,1) size(EEG.data,2) size(EEG.data,3)]);
                        else
                            EEG.icaact = EEG.icaweights*EEG.icasphere*EEG.data;
                        end
                    end
                    
                    ncomp = length(EEG.icawinv); % number of components
                    
                    % Eye blinks and saccades detection by correlation with VEOG and HEOG
                    VEOG = dataset_eog.data(1,:,:);
                    VEOG = VEOG(:);
                    HEOG = dataset_eog.data(2,:,:);
                    HEOG = HEOG(:);
                    ICs = EEG.icaact(:,:)';
                    for ic = 1:size(ICs,2)
                        corr_V(ic) = corr(ICs(:,ic),VEOG);
                        corr_H(ic) = corr(ICs(:,ic),HEOG);
                    end
                    rej_V = myFindOutliers(corr_V);
                    rej_H = myFindOutliers(corr_H);
                    Rej_ic_eog = unique([rej_V,rej_H]);     % ICs containing blinks
                    clear rej_V rej_H % free some memory
                    
                    % ICs with generics discontinutiy of spatial features
                    topography = EEG.icawinv';  % topography of the IC weigths
                    channel = EEG.chanlocs(EEG.icachansind');
                    xpos=[channel.X];ypos=[channel.Y];zpos=[channel.Z];
                    pos=[xpos',ypos',zpos'];
                    gen_disc = zeros(1,size(ICs,2)); % generic discontinuity
                    for ic = 1:ncomp
                        aux = [];
                        for el = 1:length(channel)-1
                            
                            P_el = pos(el,:); %position of current electrode
                            d = pos - repmat(P_el,length(channel),1);
                            dist = sqrt(sum((d.*d),2));
                            
                            [y,I] = sort(dist);
                            rep_ch = I(2:11); % the 10 nearest channels to el
                            weight_ch = exp(-y(2:11)); % respective weights, computed wrt distance
                            
                            aux = [aux abs(topography(ic,el)-mean(weight_ch.*topography(ic,rep_ch)'))];
                            % difference between el and the average of 10 neighbor el
                            % weighted according to weight
                        end
                        gen_disc(ic)=max(aux);
                    end
                    Rej_ic_gd = myFindOutliers(gen_disc); % ICs containing generic discontinuities
                    clear channel aux pos topography % free some memory
                    
                    % Muscle activity usually has low autocorrelation of time course
                    ncorrint =round(25/(1000/EEG.srate)); % number of samples for 25 ms lag
                    for k = 1:ncomp
                        y = EEG.icaact(k,:,:);
                        yy = xcorr(mean(y,3),ncorrint,'coeff');
                        autocorr(k) = yy(1);
                    end
                    Rej_muscle = myFindOutliers(autocorr);
                    
                    % Drop the ICs
                    IC_droplist = unique([Rej_ic_eog Rej_muscle Rej_ic_gd]); % save the number of the components dropped
                    if ~isempty(IC_droplist)
                        EEG = pop_subcomp( EEG, IC_droplist, 0); %drop components
                    end
                    EEG = eeg_checkset( EEG );
                    % ------------------------------------------------------------------------------------
                else
                    % when we do not use ICA
                    EEG = dataset;
                    IC_droplist = [];
                end
                
                clear dataset
                
                % Check if still have some artifacts after IC -- focus on each
                % channel of each trial
                EEGtmp = EEG; % local copy
                for itrial = 1:EEG.trials
                    EEGtmp.data = EEG.data(:,:,itrial);
                    EEGtmp.trials = 1;
                    EEGtmp.event = [];
                    % Mean diff value
                    [mean_diff,~]=myBiweight(diff(EEGtmp.data,[],2));
                    % Variance of the channels
                    [~,chan_std]=myBiweight(EEGtmp.data);
                    
                    % Find the outlier channels
                    Rej_mean_diff = myFindOutliers(mean_diff);
                    Rej_chan_std = myFindOutliers(chan_std);
                    %                 Rej_h_fd = myFindOutliers(h_fd);
                    % interpolation of bad epoch channel
                    %                 bad_ch_trial = unique([Rej_mean_diff Rej_chan_std]); % will be removed later
                    bad_ch_trial = intersect(Rej_mean_diff, Rej_chan_std); % 18/07/2017 %use matlab fct intersect as unique is a eeg-lab function
                    if ~isempty(bad_ch_trial)
                        EEGtmp = eeg_interp(EEGtmp, bad_ch_trial, 'spherical');
                    end
                    % add the data to the dataset
                    EEG.data(:,:,itrial) = EEGtmp.data;
                    
                    % count the number of bad channel in each trial
                    trial_bad_ch = trial_bad_ch + length(bad_ch_trial);
                    
                end
                clear EEGtmp % clear some memory
                
                %-------------------------------------------------------------------------------
                
                % remove baseline in epoched data if prest is one
                if doRemPrestBaseline == 1
                    EEG = MyRmBaseline(EEG,prestim);
                end
                
                % reference to the average
                EEG = pop_reref( EEG, []);
                EEG = eeg_checkset( EEG );
                
                %-------------------------------------------------------------------------------
                % Add epochs to corresponding conditions
                % Condition 1 - Long SOA
                l_soa_ind = find(data_dv.context_no == 1);
                try EEG_l_soatmp = pop_select(EEG,'trial',l_soa_ind); catch end;
                % Condition 2 - Short SOA
                s_soa_ind = find(data_dv.context_no == 2);
                try EEG_s_soatmp = pop_select(EEG,'trial',s_soa_ind); catch end;
                % Condition 3 - Mask only
                mask_ind = find(data_dv.context_no == 3);
                try EEG_masktmp = pop_select(EEG,'trial',mask_ind); catch end;
                % Condition 4 - Vernier only
                vernier_ind = find(data_dv.context_no == 4);
                try EEG_verniertmp = pop_select(EEG,'trial',vernier_ind); catch end;
                
                if iRun==1;
                    EEG_l_soa = EEG_l_soatmp;
                    EEG_s_soa = EEG_s_soatmp;
                    EEG_mask = EEG_masktmp;
                    EEG_vernier = EEG_verniertmp;
                    behavioral_data = data_dv;
                else
                    try EEG_l_soa = pop_mergeset(EEG_l_soa,EEG_l_soatmp); catch end;
                    try EEG_s_soa = pop_mergeset(EEG_s_soa,EEG_s_soatmp); catch end;
                    try EEG_mask = pop_mergeset(EEG_mask,EEG_masktmp); catch end;
                    try EEG_vernier = pop_mergeset(EEG_vernier,EEG_verniertmp); catch end;
                    behavioral_data = [behavioral_data, data_dv];
                end
                
                
            end
            
            cd(strcat(SaveDir,'\',char(Group{iGroup})));
            
            save(strcat(char(Subject_pool{iSubject}),'_preprocess.mat'),'EEG_l_soa','EEG_s_soa',...
                'EEG_mask','EEG_vernier','bad_trials','bad_channels','invalid_trials',...
                'blink_trials', 'trial_bad_ch','behavioral_data'); %'error_trials'
            
       catch
            cd(strcat(SaveDir,'\',char(Group{iGroup})));
            save(strcat(char(Subject_pool{iSubject}),'_error.mat'), 'SubjDir');
       end
        
    end
    
    
end



