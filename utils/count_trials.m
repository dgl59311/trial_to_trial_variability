function [ntrials] = count_trials(datadir, name, beh)
    filedir = fullfile(datadir, name);
    eeg = load(filedir);
	condition = fieldnames(eeg);
    eeg = eeg.(condition{1});
    if ~isempty(beh) && ~contains(name, 'mask')
        behdata = load(fullfile(beh, [name '_beh.mat'])).condhits;
        eeg = pop_select(eeg, 'trial', find(behdata)); 
    end    
    ntrials = eeg.trials;
end