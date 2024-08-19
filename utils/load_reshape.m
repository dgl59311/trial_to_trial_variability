function [eeg, ntrials, lentrial] = load_reshape(datadir, name, mintrials, beh)
    filedir = fullfile(datadir, name);
    eeg = load(filedir);
    condition = fieldnames(eeg);
    eeg = eeg.(condition{1});
    if ~isempty(beh) && ~contains(name, 'mask')
        behdata = load(fullfile(beh, [name '_beh.mat'])).condhits;
        eeg = pop_select(eeg, 'trial', find(behdata)); 
    end    
    ntrials = eeg.trials;
    lentrial = eeg.pnts;
    if ~isempty(mintrials)
        eeg = pop_select(eeg, 'trial', 1:mintrials);
    end
    eeg = pop_reref(eeg, []);
    eeg = eeg_epoch2continuous(eeg);
    eeg = eeg.data;
end