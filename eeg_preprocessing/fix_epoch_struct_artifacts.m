%12-mar-2018    13:00
%Samran         Auckland

% function to fix the epoch struct and count of artifacts after manual
% marking/unmarking of the epochs as artifact. Initially, the EEG.reject
% and EEG.EVENTLIST were fine, but ERPLAB uses EEG.epoch struct's flag for
% counting/showing/printing artifacts table on the command window/in the
% file, which added the manual artifacts, but didn't fixed the count of
% other artifacts, although the total accepted and rejected epochs count
% was fine. Need for fix was required to have all epochs in the dataset,
% and ERPLAB can use the clean, noise and all epochs for plotting and
% analysis.
% This function will not be required any further after the healthy SEPs
% data, because mark_epochs_gui is edited to have all the structs and
% counts correctly assigned.

%% main_processing code to call this function
% this was the code in the main_processing in the for loop which called
% this function, which loaded the marked epochs file and edited the
% structs, and resaved it overwriting the original file, alongwith the txt
% files containing EVENTLIST info and epochs/artifacts stat

% % % % % % %% epoch struct fix (temporary)

% % % % % % if (fix_epoch_struct)
% % % % % %     waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
% % % % % %         [num2str(i), '/', num2str(length(folders)), '-', 'Fixing Epoch Struct and Artifacts Info']);
% % % % % %     
% % % % % %     filepath2 = [filepath, filesep, 'Marked Epochs - ', filter_type,];
% % % % % %     
% % % % % %     tmp = {ALLEEG.setname};
% % % % % %     set_idx = find(strcmp(tmp, [filename, '-Epoched-ARTmanual']));
% % % % % %     if (isempty(set_idx))
% % % % % %         EEG = pop_loadset('filename',[filename, '-Epoched-ARTmanual.set'],'filepath',filepath2);
% % % % % %         [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
% % % % % %         EEG = eeg_checkset( EEG );
% % % % % %     else
% % % % % %         EEG = pop_newset(ALLEEG, EEG, CURRENTSET, 'retrieve', set_idx);
% % % % % %         EEG = eeg_checkset( EEG );
% % % % % %     end
% % % % % %     
% % % % % %     [ALLEEG, EEG, CURRENTSET] = fix_epoch_struct_artifacts (ALLEEG, EEG, CURRENTSET, filter_type);
% % % % % % end

%%

function [ALLEEG, EEG, CURRENTSET] = fix_epoch_struct_artifacts (ALLEEG, EEG, CURRENTSET, filter_type)

% get subject number
filename = strsplit(EEG.filename, '-');
filename = filename{1};

for i = 1:length(EEG.etc.EVENTLIST.modified)
    EEG.EVENTLIST.eventinfo(i).flag = EEG.etc.EVENTLIST.modified(i);
    EEG.epoch(i).eventflag = EEG.etc.EVENTLIST.modified(i);     % added 06-march-2018. This doesn't update when erplab2eeglab is used to syncronize, so manually overwrite this.
end

EEG.filepath = [EEG.filepath filesep];

% transfer EVENTLIST to REJECT using erplab.
EEG = pop_syncroartifacts(EEG, 'Direction', 'erplab2eeglab');

% save the event list txt for modified markers
EEG = pop_exporteegeventlist(EEG , ...
    'Filename', [EEG.filepath, filename, '-eventList-ARTmanual.txt']);


% save stats in mat file
main_folder = strsplit(EEG.filepath, filesep);
main_folder = main_folder (1:end-3);
main_folder = strjoin (main_folder, filesep);

[EEG, tprej, acce, rej, histoflags] = pop_summary_AR_eeg_detection(EEG, [EEG.filepath, filesep, filename, '-stats-ARTmanual.txt']);   %prints table to file
[~, ~, ~, ~, ~] = pop_summary_AR_eeg_detection(EEG, '');   %prints table on command line  % added 06-march-2018 since separate_manually_marked_artifacts script is commented below

load ([main_folder, filesep, 'epoch_artifacts_stats_', filter_type, '.mat']);

if (~exist('stats_ARTmanual', 'var'))
    stats_ARTmanual = [];
end

%     ALLEEG = handles.ALLEEG;                  % commented 06-march-2018
%     EEG = handles.EEG;                        % commented 06-march-2018
%     separate_manually_marked_artifacts;       % commented 06-march-2018
number_rejected_trials = length(find(EEG.etc.EVENTLIST.modified ~= 0));

stats_ARTmanual(length(stats_ARTmanual)+1).subject = filename;
stats_ARTmanual(length(stats_ARTmanual)).accepted = length(EEG.etc.EVENTLIST.modified) - number_rejected_trials;
stats_ARTmanual(length(stats_ARTmanual)).rejected = number_rejected_trials;  %total epochs rejected after visual inspection
stats_ARTmanual(length(stats_ARTmanual)).percent_rejected = 100*number_rejected_trials/length(EEG.etc.EVENTLIST.modified);
stats_ARTmanual(length(stats_ARTmanual)).F8_others = histoflags(8);

save ([main_folder, filesep, 'epoch_artifacts_stats_', filter_type, '.mat'], 'stats_ARTmanual', '-append');

% save EEG dataset. NOTE: it overwrites the dataset.
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = pop_saveset(EEG, 'savemode','resave');
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

end