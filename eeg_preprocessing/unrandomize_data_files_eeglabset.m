%19-mar-2018    16:11
%Samran         Auckland

%24-sep-2018    10:58
%Samran         Aalborg
% modified to unrandomize the EEGLAB dataformat.

% unrandomize data
% this function is to be used after analysis of EEGLAB datasets is
% completed, and unblinding is required. 
% this function loads EEGLAB datasets from filter folders, FIR and/or IIR,
% and saves them with changing EEGLAB dataset name and file name, with
% hierarchy: 
% Subject number/eeg_type/intervention/intervention-eeg_type.set
% If there are multiple sessions, pre-seps2, post-seps3, it keeps the last
% version only, with the folder and filename pre-seps.

function [] = unrandomize_data_files_eeglabset()

%% code and data folders
code_folder = cd(['..', filesep, 'Data Randomized-rest']);

data_folder = cd (code_folder);
data_folder = [data_folder, filesep];


save_folder = [data_folder, 'EEGLAB Datasets', filesep];

%% filter
filter_type = {'FIR'}; % {'FIR', 'IIR'};        % filter type: FIR or IIR. Put in cell array {'FIR', 'IIR'} if need to do for both. However, give only one when converting EEGLAB datasets to ERPLAB datasets.
filter_bandpass.HP_cutoff = 0.5;         % HP cutoff frequency in Hz: needed for saving/loading folder info
filter_bandpass.LP_cutoff = 1000;       % LP cutoff frequency in Hz: needed for saving/loading folder info

%% No ICA or with ICA
isICA = {'-ICA-'};       % give {'-NoICA-', '-ICA-'} to do both 'only epoch cleaned with No ICA' and 'epoch cleaned+ICA'. This is overwritten when doing EEGLAB dataset to ERPLAB dataset conversion.


%% load data_secret.mat
load ([data_folder, 'data_info', filesep, 'data_secret.mat']);

% remove noisey alzheimers rest
tmp = [data_info.subject];  % alzheimers, remove subject 7 as already analyzed as stroe subject 4.
tmp = find(tmp == 7 | tmp == 8 | tmp == 12 | tmp == 14);
data_info(tmp) = [];


%%
w = waitbar (0, '', 'Name', 'UnRandomizing files');
waitbar_steps = length(data_info) * length(filter_type) * length(isICA);
waitbar_step = 0;
waitbar_resolution = 1/1;       %1 divided by the number of times waitbar is to be updated in one loop (inner)

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for di = 1:length(data_info)
    random_number = data_info(di).random_number;
%     subject = str2num(data_info(di).subject);
    subject = data_info(di).subject;
    setname = data_info(di).setname;
    
    setname = strsplit (setname, '-');
    
    intervention = setname{2};        % Aerobics, Ctrl, Chiro
    session = setname{3};             % pre, post
    eeg_type = setname{4};            % seps, rest, seps2
    
    % if multiple sessions, keep only the last one (the one with highest
    % integer). the data_info contains setname in sequence, therefore,
    % it'll be pre-seps, then pre-seps2, so on. They are in order. So, the
    % following codes works.
    % likely the highest integer one is the better session, however, still
    % verify in the preprocessing, and then modify here if required if any
    % other session for a particular subject is better than the highest
    % integer one.
    eeg_type_letters = find (isletter(eeg_type));
    eeg_type = eeg_type (eeg_type_letters(1):eeg_type_letters(end));
    
    for f = 1:length(filter_type)
        
        filepath = ['Subject ', num2str(subject), filesep, ...
            eeg_type, filesep, intervention, filesep, ...
            session, '-', eeg_type];
        
        save_filepath = [save_folder, 'sessions ', cell2mat(filter_type(f)), '(', num2str(filter_bandpass.HP_cutoff), '-', num2str(filter_bandpass.LP_cutoff), ')', filesep, filepath];
        
        if (~exist(save_filepath))
            mkdir (save_filepath);
        end
        
        for I = 1:length(isICA)
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                    [num2str(di), '/', num2str(length(data_info)), '-', 'Unrandomizing EEGLAB dataset - ', cell2mat(filter_type(f)), cell2mat(isICA(I))]);
                
            if (strcmp(cell2mat(isICA(I)), '-ICA-'))
                load_filepath = [data_folder, 'Subject ', num2str(subject), filesep, num2str(random_number), filesep, 'ICA-cleaned - ', cell2mat(filter_type(f)), '(', num2str(filter_bandpass.HP_cutoff), '-', num2str(filter_bandpass.LP_cutoff), ')'];
                
                filename = [num2str(random_number), '-', 'ICA-cleaned', '.set'];
                
            elseif ((strcmp(cell2mat(isICA(I)), '-NoICA-')))
                load_filepath = [data_folder, 'Subject ', num2str(subject), filesepnum2str(random_number), filesep, 'Epochs-cleaned - ', cell2mat(filter_type(f)), '(', num2str(filter_bandpass.HP_cutoff), '-', num2str(filter_bandpass.LP_cutoff), ')'];
                
                filename = [num2str(random_number), '-', 'Epochs-cleaned', '.set'];
            end
            
            EEG = pop_loadset('filename', filename, 'filepath', load_filepath);
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            EEG = eeg_checkset( EEG );
            EEG.data = double(EEG.data);
            
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'setname',[session, '-', eeg_type],'savenew',[save_filepath, filesep, session, '-', eeg_type, '.set'],'gui','off');
            EEG.data = double(EEG.data);
            
            % reset the variables for next iteration/file
            STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        end
    end
end

%close EEGLAB gui
h = findobj('tag','EEGLAB');
close (h);

close (w);
% clc

clear

end