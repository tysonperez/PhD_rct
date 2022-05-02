%16-oct-2017    23:41
%Samran         Auckland

% verify number of epochs and their inter-duration
% manually go through the saved file struct, enter the false epoch number
% and see the neighbouring epoch timings in the EEGLAB edit->event values, 
% and on channel scroll plots to
% verify if the trigger(s) are correct. It can be that only one false_epoch
% is detected, e.g. if 20 triggers were given, they were fine, but stopped
% and triggering started again after sometime. In this case, the first 20
% triggers may have to be deleted (e.g. if total triggers are more than 
% 1000). 
% For deletion, modify the false_epoch
% values in the struct after verification manually, e.g. the triggers wrong
% are 1,2 and 1003, but the false_epochs will have the values 1,2 and 1002,
% as this is the way they are calculated by this script. Therefore, modify
% them manually to [1,2,1003] to be given as input to EEGLAB scripts to
% automatically delete them before preprocessing, and overwrite the
% epoch_info.mat file.

function [number_of_epochs, false_epochs] = check_false_triggers (EEG)

stimulation_frequency = 2.3;        % SEPs stimulation frequency in Hz
percentage_error = 25;              % percentage of stimulation frequency to consider in stimulation range

min_freq = stimulation_frequency - (stimulation_frequency * percentage_error)/100;
max_freq = stimulation_frequency + (stimulation_frequency * percentage_error)/100;

min_diff = round (EEG.srate/max_freq);      % minimum samples difference
max_diff = round (EEG.srate/min_freq);      % maximum samples difference

epochs = cell2mat({EEG.event.latency});

number_of_epochs = length(epochs);

epochs_diff = diff (epochs);        % difference between samples of consecutive triggers

false_epochs = find ((epochs_diff < min_diff) | (epochs_diff > max_diff));

end