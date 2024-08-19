%% Get data for each condition from preprocessed files
clc;
clear;
CurrDir = cd();

%%
% get data directories
DataDir         = uigetdir([], 'path where preprocessed data are stored');
ResultsFolder   = uigetdir([], 'path to store data by condition'); % use local folder
DemogDir        = uigetdir([], 'path where demographical data are stored');

%% store individual conditions in a folder and reject high amplitude trials 
groups = {'patients', 'controls'};
conditions = {'EEG_l_soa', 'EEG_s_soa', 'EEG_mask', 'EEG_vernier'};
contextn = [1, 2, 3, 4];

%%
ii = 1;
totalIterations = 960;
h = waitbar(0, 'Processing...');
for i = 1:length(groups) 
    Participants = dir(fullfile(DataDir, groups{i}, '*.mat*'));
    for j = 1:length(Participants)
        id_p = Participants(j).name;
        data = load(fullfile(DataDir, groups{i}, id_p));
        for k = 1:length(conditions)
            f_name = strrep(id_p, 'preprocess', conditions{k});
            f_name = strrep(f_name, 'EEG', groups{i}(1:3));
            cond_data = data.(conditions{k});
            hitsdata = data.behavioral_data;
            strsize = size(hitsdata, 2);
            condhits = [];
            for jj = 1:strsize
                tmpbeh = hitsdata(jj);
                tmpcond = find(tmpbeh.context_no == k);
                tmphits = tmpbeh.hits(tmpcond);
                condhits = [condhits tmphits'];
            end
            % reject trials where the data exceeded 90 and -90 uV
            before_90 = cond_data.trials;
            [cond_data, Indexes] = pop_eegthresh(cond_data, 1, 1:cond_data.nbchan,...
                -90, 90, cond_data.xmin, cond_data.xmax, 0, 1); 
            try
                condhits(Indexes) = [];
                if k ~= 3
                    save(fullfile(CurrDir, 'behavior_eeg', [f_name '_beh.mat']), 'condhits');
                end
            catch
                condhits = [];
                if k ~= 3
                    save(fullfile(CurrDir, 'behavior_eeg', [f_name '_ERROR_beh.mat']), 'condhits');
                    % behavior and brain data did not match
                end
            end
            
            save(fullfile(ResultsFolder, f_name), 'cond_data');
            
            % save info preprocessing
            info_trials{ii, 1} = f_name;
            info_trials{ii, 2} = conditions{k};
            info_trials{ii, 3} = before_90;
            info_trials{ii, 4} = cond_data.trials;
            info_trials{ii, 5} = length(condhits);
            ii = ii + 1;
            disp(ii)
            clear cond_data
            clear condhits
            
            % Update waitbar
            waitbar(ii / totalIterations, h, sprintf('Processing... %d%%', round(ii / totalIterations * 100)));
        end
    end
end
%%
% write csv table
main_results = array2table(info_trials);
main_results.Properties.VariableNames = {'ID', 'condition', ...
    'init trials', 'available trials', 'beh trials'};
writetable(main_results,...
    fullfile(DemogDir, 'available_EEG_trials.csv'));

%%