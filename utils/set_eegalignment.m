function [EEG] = set_eegalignment(data, channels, times, ntrials, srate)
    EEG.setname     = 'For alignment';
    EEG.data        = reshape(data, [1 size(data)]);
    EEG.icaact      = [];
    EEG.icawinv     = [];
    EEG.icasphere   = [];
    EEG.icaweights  = [];
    EEG.pnts        = size(data, 1);
    EEG.srate       = srate; 
    EEG.trials      = ntrials;
    EEG.nbchan      = 1;
    EEG.comments    = ' ';
    EEG.chanlocs    = channels;
    EEG.filename    = 'x';
    EEG.filepath    = '';
    EEG.times       = times;
    EEG.xmin        = times(1);
    EEG.xmax        = times(end);
    
end