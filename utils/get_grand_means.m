function [grand_average_struct] = get_grand_means(data, group, condition, behdir)
    
    [condition_data] = find_in_cell(data, condition);
    condition_data = condition_data{1};
    
    for i = 1:length(condition_data)
        eeg_file = load(fullfile(condition_data(i).folder, condition_data(i).name));
        condition_id = fieldnames(eeg_file);
        tmpdata = eeg_file.(condition_id{1});
        if ~isempty(behdir)
            try
                hits = load(fullfile(behdir, [condition_data(i).name '_beh.mat'])).condhits;
                tmpdata = pop_select(tmpdata, 'trial', find(hits));
            catch
                disp(['eeg and behavior files do not match: Check subject N ' num2str(i)])
            end
        end
        average_trials = mean(tmpdata.data, 3);
        % We convert to eeglab to remove baseline and re-reference
        eeglab_struct = set_eeglabstruct(eeg_file.(condition_id{1}), average_trials);
        concatenate_data(:, :, i) = eeglab_struct.data;
    end
    
    grand_average_struct = set_eeglabstruct(eeg_file.(condition_id{1}), mean(concatenate_data, 3));
    grand_average_struct.setname = ['Grand average ' group{1} condition{1}];
end