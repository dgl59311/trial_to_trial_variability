function [times] = get_timesvec(filedir)
    eeg = load(filedir);
    condition = fieldnames(eeg);
    eeg = eeg.(condition{1});
    times = eeg.times;
end