% script to get the number of manual selected artifact epochs. Otherwise,
% the EEG.etc.EVENTLIST.reject_difference, or
% EEG.etc.EVENTLIST.modified (bit 8) can be used to separate these.
% The following way was hit/trial (and worked) to get the number of epochs
% manually selected as artifacts for the count/stats etc.

setname = EEG.setname;
EEG  = pop_resetrej(EEG , 'ResetArtifactFields', 'off', 'ArtifactFlag',  1:7);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname', [setname, '-manual-selected-only'], 'gui','off');

EEG = pop_syncroartifacts(EEG, 'Direction', 'erplab2eeglab');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'overwrite', 'on', 'gui','off');

[EEG, tprej, acce, rej, histoflags ] = pop_summary_AR_eeg_detection(EEG, '');   %prints table on command line

clear setname;

eeglab redraw;