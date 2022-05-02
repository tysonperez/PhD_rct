
function [] = main_preprocessing()

addpath ('E:\Folders\PhD\ISAD\ISAD_Curry_data\ISAD_EEGs\Preprocessing\eeg_preprocessing\eeglab14_1_1b\eeglab14_1_1b');
eeglab

%% code and data folders
  code_folder = cd(['E:\Folders\PhD\ISAD\ISAD_Curry_data\ISAD_EEGs\Preprocessing',filesep,'Healthy']);

data_folder = cd (code_folder);
data_folder = [data_folder, filesep];

%% eeg type
% eeg_type = {'rest', 'seps'};
eeg_type = {'rest'};

%% steps to perform: 1 & 3 true, all others false > 4 true, all others false > 5 true, all other false > etc

process_raw_file = true;               % process the raw file. this is done in combination with either checking_false_triggers or fix_channels_and_triggers. Needs to be done before running PREP and further pipeline processes.
checking_false_triggers = false;        % ALWAYS FALSE
fix_channels_and_triggers = true;      % after identification and correction of faulty triggers in the epoch_info, this will remove redundant channels, add channel info, remove faulty triggers and remove redundant data in the beginning and end. saves 0th file.
run_PREP = false;                       % run PREP pipeline
generate_event_list = false;            % generate EventList using ERPLAB
mark_epochs_rejection = false;          % mark epochs automatically using ERPLAB
run_ICA = false;                        % run ICA
ICA_with_PCA = false;                   % ALWAYS FALSE
apply_IC_weights_rejection = false;     % true: filter, epoch and apply ICA weights to save dataset ready for manual marking of ICs for rejection. calculate other parameters to help in marking ICs.
reject_ICA_components = false;          % true: reject marked ICs (need to mark ICs manually and verify how data looks, save it and then run codes with this as true)


%% Post-ICA Filter (FIR then IIR)
% This filter is used after running ICA (i.e after mark_epochs_rejection and run_ICA). 
% The IC weights are applied to FIR filtered (0.5-100 Hz, trans band 1 Hz, order 3624) data for marking ICs. 
% Afterwards, the IIR filtered (Butterworth 0.01-100 Hz, order 1) is applied when ICs are rejected giving the cleaned dataset on which the analysis are done.
% These filters are different from the filter in the next section applied for marking epochs and ICA calculations. 
% The filter_type (i.e., FIR) set here for applying ICA weights/matrix is
% also used for marking epochs & ICA but the cutoff frequencies, transition band, order are different.
% Finally, when reject_ICA_components is run, will automatically switch to from FIR to IIR in the code but remember to change filter_order to 1 and bandpass to 0.01-100 Hz

if (strcmp(eeg_type, 'seps'))
    filter_type = 'IIR';                   
    filter_order = 2;                       
    filter_bandpass.HP_cutoff = 0.5;        
    filter_bandpass.LP_cutoff = 1000;       
elseif (strcmp (eeg_type, 'rest'))
    filter_type = 'FIR';             % keep FIR here for loading and saving files, but for the final step the code will filter IIR as the if/else check condition is removed                  
    filter_order = 3624;             % apply_IC_weights_rejection = 3624; reject_ICA_components =1; for determining FIR order (ie. 3624) when applying ICA weights ran pop_firwsord() function from Matlab command window (samp rate 1000 Hz, kaiser 5.65, max ripple 0.001, trans band 1 Hz)          
    filter_bandpass.HP_cutoff = 0.5; % apply_IC_weights_rejection = 0.5; reject_ICA_components = 0.01         
    filter_bandpass.LP_cutoff = 100;            
end

%% Marking Epochs & ICA filter (FIR)

ICA_downsample_rate = 500; % keep in mind that total amount data should be 20 to 30 * # of channels^2. Apply the ICA weights calculated with this sampling rate to the original sampling rate data.
ICA_temp_folder = 'E:\Folders\PhD\ISAD\ISAD_Curry_data\ISAD_EEGs\Preprocessing\eeg_preprocessing\temp_AMICA';      %temporary folder for AMICA to load and store files, as it doesn't work with spaces in the folder, file names. make this folder manually before running AMICA

% filter_ICA is used in HP (and LP if required) filtering in marking of
% epochs, and the saved HP filtered file is also used/loaded when running
% ICA to save time in filtering the data again.

if (strcmp(eeg_type, 'seps'))
    filter_ICA_order = 2;           
    filter_ICA.HP_cutoff = 1;      
    filter_ICA.LP_cutoff = [];      
elseif (strcmp(eeg_type, 'rest'))
    filter_ICA_order = 2416;        % pop_firwsord() = FIR order for kaiser 5.65, max ripple 0.001, transition band 1.5Hz at 1000 Hz sampling rate.
    filter_ICA.HP_cutoff = 1;       % HP cutoff frequency in Hz
    filter_ICA.LP_cutoff = [];      % LP cutoff frequency in Hz. Default is empty. But after HP filtering, if epochs are difficult to identify as artifcats, HP and LP the data (Bandpassed) to aid in visualization of epochs.
end

%% load basic files

load ([data_folder, 'data_info', filesep, 'data_share.mat']);     % file containing information shared as the data is blinded. file data_secret contains the complete information of the files to unblind.

% data_info([5, 20]) = [];                                    % subject not to be included in pipeline (5: PREPed didn't work because of too noisy data; 20: incomplete data)
% data_info([1:23,25:end]) = [];

% subject_id = cell2mat({data_info.subject});
subject_id = str2double({data_info.subject});       % for Tyson data_info character type subject_id


%% calculate number of steps for waitbar, based on the processing steps being performed
waitbar_step_counter = 0;
if (process_raw_file)
    waitbar_step_counter = waitbar_step_counter + 1;
end

if (checking_false_triggers)
    waitbar_step_counter = waitbar_step_counter + 1;
end

if (fix_channels_and_triggers)
    waitbar_step_counter = waitbar_step_counter + 2;
end

if (run_PREP)
    waitbar_step_counter = waitbar_step_counter + 1;
end

if (generate_event_list)
    waitbar_step_counter = waitbar_step_counter + 1;
end

if (mark_epochs_rejection)
    waitbar_step_counter = waitbar_step_counter + 3;
    
    if (~isempty(filter_ICA.LP_cutoff))
        waitbar_step_counter = waitbar_step_counter + 1;
    end
    
    epoch_rej_parameters = epoch_rejection_parameters (length(unique(subject_id)));
end

if (run_ICA)
    waitbar_step_counter = waitbar_step_counter + 5;
end

if (apply_IC_weights_rejection)
    waitbar_step_counter = waitbar_step_counter + 5;
end

if (reject_ICA_components)
    waitbar_step_counter = waitbar_step_counter + 7;
end

w = waitbar (0, '', 'Name', 'Preprocessing');
waitbar_steps = length (subject_id);
waitbar_step = 0;
waitbar_resolution = 1/waitbar_step_counter;       % 1 divided by the number of times waitbar is to be updated in one loop (inner)


%% start

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

if (process_raw_file)   % modify the options of EEGLAB, assuming process_raw_file is run before analysis of a new data starts, otherwise modify manually, especially change single precision to double precision
    pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 1, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 0);
end

missing_files = [];

for i = 1:length(subject_id) 
    try
        
        filepath = [data_folder, 'Subject ', data_info(i).subject, filesep, data_info(i).session];
        
        filename = data_info(i).session;
        
        %% Raw file EEGLAB
        
        if (process_raw_file)
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Loading File']);
            
            %% load file
            EEG = pop_loadset('filename',[filename, '-Raw.set'],'filepath',[filepath, filesep, 'Raw']);
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            EEG = eeg_checkset( EEG );     
            
             %% Remove redundant channels, add channel info, and remove redundant data in the start and end
            if (fix_channels_and_triggers)
                waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                    [num2str(i), '/', num2str(length(subject_id)), '-', 'Removing redundant ch & adding ch info']);
                
                %if (isempty(EEG.chanlocs))          %if the EEGLAB file is missing channel info
                    %load ('MNI_60ch_locs.mat');
                    
                    %EEG.chanlocs = chanlocs;
                    %EEG.etc.missing_channel_locations = true;
                %end
                
                if (EEG.nbchan == 69)
                    EEG = pop_select( EEG,'nochannel',{'VEOG' 'HEOG' 'EKG' 'EMG' 'GSR' 'F11' 'FT11' 'F12' 'FT12'});
                else
                    EEG = pop_select( EEG,'nochannel',{'VEOG' 'HEOG' 'EKG' 'EMG' 'GSR' 'F11' 'FT11' 'F12' 'FT12','SAW'});
                end
                [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[filename, '-60ch'],'gui','off');
                EEG = eeg_checkset( EEG );
                
                %%% add MNI coordinate file here (not BESA!)
                %% to get the code EEG=pop_chanedit() below follow these steps...
                %% EEGLAB GUI select file > load existing dataset (raw folder .set file) > edit > select data > remove ch without coordinates > Ok > edit > channel locations > change channel names if needed (e.g. CB1 to I1) > look up locs > MNI > optimize head center > Ok > type eegh on matlab command line  
                EEG=pop_chanedit(EEG, 'changefield',{59 'labels' 'I1'},'changefield',{60 'labels' 'I2'},'settype',{'1:60' 'EEG'},'lookup','E:\\Folders\\PhD\\ISAD\\ISAD_Curry_data\\ISAD_EEGs\\Preprocessing\\eeg_preprocessing\\eeglab14_1_1b\\eeglab14_1_1b\\plugins\\dipfit2.3\\standard_BEM\\elec\\standard_1005.elc','eval','chans = pop_chancenter( chans, [],[]);');

                [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
                EEG = eeg_checkset( EEG );
                
%               EEG = pop_select( EEG,'nochannel',{'M1' 'M2'})
%               [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',[filename, '-58ch'],'gui','off');
%                EEG = eeg_checkset( EEG );
                
                
                
                %% remove redudant data in the beginning and end
                if (strcmp (eeg_type, 'seps'))
                    %%%% remove redudant data in the beginning and end (event-related data)
                    % if there is more than 30 sec of data before the 1st trigger or
                    % after the last trigger, it is removed to reduce the size of
                    % datasets, reduce pre-processing time, and e.g. for PREP, not
                    % making it wrong decisions for interpolation based on the redudant
                    % data.
                    
                    event_time = cell2mat({EEG.event.latency})/EEG.srate;   %in seconds
                    data_start = event_time(1) > 30;
                    data_end = EEG.xmax > event_time(end) + 30;
                    
                    if (data_start)
                        data_start = event_time(1) - 30;
                    else
                        data_start = EEG.xmin;
                    end
                    
                    if (data_end)
                        data_end = event_time(end) + 30;
                    else
                        data_end = EEG.xmax;
                    end
                elseif (strcmp (eeg_type, 'rest'))
                    %%%% remove redudant data in the beginning and end (resting-state data)
                    % take the middle 600sec of the data if data is more than 10 min
                    % long. Done so redundant data is removed to reduce the size of
                    % datasets, reduce pre-processing time, and e.g. for PREP, not
                    % making it wrong decisions for interpolation based on the redundant
                    % data.
                    
                    % if <= 600 sec data, include approx 30 secs on both
                    % ends to avoid filtering artefacts
                    
                    if (EEG.xmax > 600)
                        data_start = (EEG.xmax - 600)/2;
                        data_end = data_start + 600;
                    else
                        data_start = EEG.xmin + 30;
                        data_end = EEG.xmax - 30;
                    end
                end
                
                EEG = pop_select( EEG,'time',[data_start, data_end] );
                
                %% save 0 file, with channel data
                
                waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                    [num2str(i), '/', num2str(length(subject_id)), '-', 'Saving 0th file']);
                
                save_filepath = [filepath, filesep, '0'];
                if (~exist(save_filepath, 'dir'))
                    mkdir (save_filepath);
                end
                
                [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'setname',[filename, '-0'],'savenew',[save_filepath, filesep, filename, '-0.set'],'gui','off');
                EEG.data = double(EEG.data);    % pop_newset with save, saves in single format as double precision is not required in saving.
                [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
                EEG = eeg_checkset( EEG );
                
                ALLEEG = pop_delset( ALLEEG, [1:length(ALLEEG)] );   %save memory
            end
        end
        
        
        %% PREP pipeline
        if (run_PREP)
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'PREP', '-started:', char(datetime)]);
            
            % check if 0th file is loaded, as PREP will process that, e.g. if
            % fix_channels_and_triggers was performed earlier, and no longer
            % required to repeat and is set to false, then 0th file will not be
            % created and not present in current datasets.
            set_idx = [];
            if (~isempty(ALLEEG))
                tmp = {ALLEEG.setname};
                set_idx = find(strcmp(tmp, [filename, '-0']));
            end
            if (isempty(set_idx))
                load_filepath = [filepath, filesep, '0'];
                EEG = pop_loadset('filename',[filename, '-0.set'],'filepath',load_filepath);
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                EEG = eeg_checkset( EEG );
            else
                EEG = pop_newset(ALLEEG, EEG, CURRENTSET, 'retrieve', set_idx);
                EEG = eeg_checkset( EEG );
            end
            
            save_filepath = [filepath, filesep, 'PREPed'];
            if (~exist(save_filepath, 'dir'))
                mkdir (save_filepath);
            end
            
            EEG.data = double(EEG.data);
            
            EEG = pop_prepPipeline(EEG, struct('lineFrequencies', [50:50:EEG.srate/2], ...
                'reportMode', 'normal', ...
                'sessionFilePath', [save_filepath, filesep, 'Report.pdf'], ...
                'summaryFilePath', [save_filepath, filesep, 'Summary.html'], ...
                'consoleFID', 1));
            
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'setname',[filename, '-PREPed'],'savenew',[save_filepath, filesep, filename, '-PREPed.set'],'gui','off');
            EEG.data = double(EEG.data);    % pop_newset with save, saves in single format as double precision is not required in saving.
        end
        
        
        %% Generate event list using ERPLAB
        % get events from EEGLAB dataset and put them in BINs. Updates the
        % EEG.event structure as well. assumes events are imported previously
        % using EEGLAB and trigger channel (any of the raw, 0th, PREP file will
        % work). For present study/case, PREP file is expected and saved file
        % name contains 'PREPed' as well.
        
        if (generate_event_list)
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Generating Event List']);
            
            % check if PREP file is loaded, e.g. if PREP pipeline was run
            % earlier and no longer required to repeat and is set to false,
            % then PREP file will not be created and not present in current
            % datasets.
            set_idx = [];
            if (~isempty(ALLEEG))
                tmp = {ALLEEG.setname};
                set_idx = find(strcmp(tmp, [filename, '-PREPed']));
            end
            if (isempty(set_idx))
                load_filepath = [filepath, filesep, 'PREPed'];
                EEG = pop_loadset('filename',[filename, '-PREPed.set'],'filepath',load_filepath);
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                EEG = eeg_checkset( EEG );
            else
                EEG = pop_newset(ALLEEG, EEG, CURRENTSET, 'retrieve', set_idx);
                EEG = eeg_checkset( EEG );
            end
            
            % add events in the dataset for epoching/cleaning later for resting state EEG
            if (strcmp (eeg_type, 'rest'))
                EEG_length = 600;
                [EEG, EEG.etc.epoch_length] = add_events (EEG, EEG_length);
                EEG = eeg_checkset( EEG );
            end
            
            save_filepath = [filepath, filesep, 'eventListed'];
            if (~exist(save_filepath, 'dir'))
                mkdir (save_filepath);
            end
            
            % check elist-equations file for correctness (e.g. seps and
            % rest has different event-markers)
            EEG  = pop_editeventlist(EEG, ...
                'AlphanumericCleaning', 'on', ...
                'BoundaryNumeric', { -99}, ...
                'BoundaryString', { 'boundary' }, ...
                'ExportEL', [save_filepath, filesep, filename, '-eventList.txt'], ...
                'List', [data_folder, 'data_info', filesep, 'elist-equations-', cell2mat(eeg_type), '.txt'], ...
                'SendEL2', 'EEG&Text', ...
                'UpdateEEG', 'codelabel', ...
                'Warning', 'off' );
            
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'setname',[filename, '-eventListed'],'savenew',[save_filepath, filesep, filename, '-eventListed.set'],'gui','off');
            EEG.data = double(EEG.data);    % pop_newset with save, saves in single format as double precision is not required in saving.
            EEG = eeg_checkset( EEG );
        end
        
        
        %% Mark epochs for rejection
        
        if (mark_epochs_rejection)
            if (~exist([filepath, filesep, 'Marked Epochs - ', filter_type], 'dir'))
                mkdir ([filepath, filesep, 'Marked Epochs - ', filter_type]);
            end
            
            if (~exist('stats_ARTautomatic', 'var'))
                stats_ARTautomatic = [];
                
                if (exist ([data_folder, filesep, 'data_info', filesep, 'epoch_artifacts_stats_', filter_type, '.mat'], 'file'))
                    load ([data_folder, filesep, 'data_info', filesep, 'epoch_artifacts_stats_', filter_type, '.mat']);
                end
            end
            
            %% filter data
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Highpass Filtering']);
            
            % check if file "with event list" is loaded, e.g. if PREP
            % pipeline was run earlier, and event list created,
            % and no longer required to repeat and is set to false,
            % then PREP file "with event list" will not be created and not
            % present in current datasets.
            set_idx = [];
            if (~isempty(ALLEEG))
                tmp = {ALLEEG.setname};
                set_idx = find(strcmp(tmp, [filename, '-eventListed']));
            end
            if (isempty(set_idx))
                load_filepath = [filepath, filesep, 'eventListed'];
                EEG = pop_loadset('filename',[filename, '-eventListed.set'],'filepath',load_filepath);
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                EEG = eeg_checkset( EEG );
            else
                EEG = pop_newset(ALLEEG, EEG, CURRENTSET, 'retrieve', set_idx);
                EEG = eeg_checkset( EEG );
            end
            
            if (strcmp(filter_type, 'FIR'))
                % the order of the filter was calculated using pop_firwsord, or
                % EEGLAB->Tools->Filter the data->Windowed Sinc FIR Filter
                EEG = pop_firws(EEG, 'fcutoff', filter_ICA.HP_cutoff, 'ftype', 'highpass', 'wtype', 'kaiser', 'warg', 5.65326, 'forder', filter_ICA_order, 'minphase', 0);
            elseif (strcmp (filter_type, 'IIR'))
                [b, a] = butter (filter_ICA_order, filter_ICA.HP_cutoff/(EEG.srate/2), 'high');
                EEG.data = filtfilt (b, a, EEG.data')';
            end
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6,'setname',[filename, '-HPfiltered'],'savenew',[filepath, filesep, 'Marked Epochs - ', filter_type, filesep, filename, '-HPfiltered.set'],'gui','off');
            EEG.data = double(EEG.data);    % pop_newset with save, saves in single format as double precision is not required in saving.
            EEG = eeg_checkset( EEG );
            
            if (~isempty(filter_ICA.LP_cutoff))     % for better visualization if highpass filtered is difficult to identify/verify bad epoch visually. Lowpass filter the data and save it as bandpassed data. The data used for further analysis from this dataset is the marked epochs.
                waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                    [num2str(i), '/', num2str(length(subject_id)), '-', 'Lowpass Filtering']);
                
                if (strcmp(filter_type, 'FIR'))
                    % the order of the filter was calculated using pop_firwsord, or
                    % EEGLAB->Tools->Filter the data->Windowed Sinc FIR Filter
                    EEG = pop_firws(EEG, 'fcutoff', filter_ICA.LP_cutoff, 'ftype', 'lowpass', 'wtype', 'kaiser', 'warg', 5.65326, 'forder', filter_ICA_order, 'minphase', 0);
                elseif (strcmp (filter_type, 'IIR'))
                    [b, a] = butter (filter_ICA_order, filter_ICA.LP_cutoff/(EEG.srate/2), 'low');
                    EEG.data = filtfilt (b, a, EEG.data')';
                end
                [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6,'setname',[filename, '-BPfiltered'],'savenew',[filepath, filesep, 'Marked Epochs - ', filter_type, filesep, filename, '-BPfiltered.set'],'gui','off');
                EEG.data = double(EEG.data);    % pop_newset with save, saves in single format as double precision is not required in saving.
                EEG = eeg_checkset( EEG );
            end
            
            
            %% extract epochs
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Extracting Epochs']);
            
            if (strcmp (eeg_type, 'rest'))
                EEG = pop_epochbin(EEG, [0.0 1000*(EEG.etc.epoch_length-1/(2*EEG.srate))], 'none');
            elseif (strcmp (eeg_type, 'seps'))
                EEG = pop_epochbin(EEG, [-100.0 150.0], 'pre');     
            end
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 7,'setname',[filename, '-Epoched'],'savenew',[filepath, filesep, 'Marked Epochs - ', filter_type, filesep, filename, '-Epoched.set'],'gui','off');
            EEG.data = double(EEG.data);    % pop_newset with save, saves in single format as double precision is not required in saving.
            
            
            %% mark epochs for rejection using ERPLAB
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Marking Epochs for Rejection']);
            
            rej_parameters = epoch_rej_parameters(subject_id(i));
            % rej_parameters = epoch_rej_parameters(1);  
            % do not reject epochs here. Just mark, and after visual
            % verification of performance of automated detection of noisy
            % epochs, store the noisy epoch indices, so they can be rejected
            % whenever required, e.g. before running ICA, and before applying ICA
            % weights after ICs are computed. Don't reject, as it in future
            % all epochs may be needed, or only the marked noisy ones to be
            % plotted to see the difference of noisy vs cleaned.
            
            % the datasets created contain the information of previous type of
            % artifacts, as the file being used to detect artifacts is not the
            % raw epoched file, as the raw epoched file sequentially goes
            % through : Max Volatge, peak-to-peak Voltage, Step-Like,
            % Sample-to-Sample Difference and Flatline artifacts.
            % Stores same data in 2 files named automatic and manual. Manual
            % file when modified can be compared with the automatic one to see
            % the automated vs automated+manual artifact marking of epochs.
            
            % NOTE: modify the thresholds and windows, they are hard-coded. For
            % a study, go through with ERPLAB gui to get these values, and then
            % modify here. Modification of size of epoch (in extract epochs)
            % may also be needed.
            
            % extreme/max voltage threshold, positive and negative. Marker bits 1 and 2.
            if (strcmp (eeg_type, 'rest'))
                EEG  = pop_artextval( EEG , 'Channel',  1:EEG.nbchan, 'Flag', [1 2], 'Threshold', [-1 1]*rej_parameters.max_voltage_threshold.value, 'Twindow', [EEG.times(1) EEG.times(end)] ); 
                EEG = pop_syncroartifacts(EEG, 'Direction', 'erplab2eeglab');
            elseif (strcmp (eeg_type, 'seps'))
                EEG  = pop_artextval( EEG , 'Channel',  1:EEG.nbchan, 'Flag', [1 2], 'Threshold', [-1 1]*rej_parameters.max_voltage_threshold.value, 'Twindow', [EEG.times(1) rej_parameters.stimulus.window(1)] ); 
                EEG = pop_syncroartifacts(EEG, 'Direction', 'erplab2eeglab');
                
                EEG  = pop_artextval( EEG , 'Channel',  1:EEG.nbchan, 'Flag', [1 2], 'Threshold', [-1 1]*rej_parameters.max_voltage_threshold.value, 'Twindow', [rej_parameters.stimulus.window(2) EEG.times(end)] ); 
                EEG = pop_syncroartifacts(EEG, 'Direction', 'erplab2eeglab');
            end
            
            % peak-to-peak volage threshold in a given window. Marker bits 1 and 3.
            EEG  = pop_artmwppth( EEG , 'Channel',  1:EEG.nbchan, 'Flag', [1 3], 'Threshold',  rej_parameters.pp_threshold.value, 'Twindow', [EEG.times(1) EEG.times(end)], 'Windowsize',  rej_parameters.pp_threshold.window_size, 'Windowstep',  rej_parameters.pp_threshold.window_step );
            EEG = pop_syncroartifacts(EEG, 'Direction', 'erplab2eeglab');
            
            % step-like artifacts. Mark flag 1 and 4.
            EEG  = pop_artstep( EEG , 'Channel',  1:EEG.nbchan, 'Flag', [1 4], 'Threshold',  rej_parameters.step_like.value, 'Twindow', [EEG.times(1) EEG.times(end)], 'Windowsize',  rej_parameters.step_like.window_size, 'Windowstep',  rej_parameters.step_like.window_step );
            EEG = pop_syncroartifacts(EEG, 'Direction', 'erplab2eeglab');
            
            % sample-to-sample difference artifacts. Marker bits 1 and 5.
            if (strcmp (eeg_type, 'rest'))
                EEG  = pop_artdiff( EEG , 'Channel',  1:EEG.nbchan, 'Flag', [1 5], 'Threshold',  rej_parameters.sample_difference.value, 'Twindow', [EEG.times(1) EEG.times(end)] );
                EEG = pop_syncroartifacts(EEG, 'Direction', 'erplab2eeglab');
            elseif (strcmp (eeg_type, 'seps'))
                EEG  = pop_artdiff( EEG , 'Channel',  1:EEG.nbchan, 'Flag', [1 5], 'Threshold',  rej_parameters.sample_difference.value, 'Twindow', [EEG.times(1) rej_parameters.stimulus.window(1)] );
                EEG = pop_syncroartifacts(EEG, 'Direction', 'erplab2eeglab');
                
                EEG  = pop_artdiff( EEG , 'Channel',  1:EEG.nbchan, 'Flag', [1 5], 'Threshold',  rej_parameters.sample_difference.value, 'Twindow', [rej_parameters.stimulus.window(2) EEG.times(end)] );
                EEG = pop_syncroartifacts(EEG, 'Direction', 'erplab2eeglab');
            end
            
            % flatline artifacts. Marker bits 1 and 6.
            EEG  = pop_artflatline( EEG , 'Channel',  1:EEG.nbchan, 'Duration',   rej_parameters.flatline.duration, 'Flag', [1 6], 'Threshold', [-1 1]*rej_parameters.flatline.value, 'Twindow', [EEG.times(1) EEG.times(end)] );
            EEG = pop_syncroartifacts(EEG, 'Direction', 'erplab2eeglab');
            
            % save
            EEG.etc.EVENTLIST.original = EEG.EVENTLIST;     %store EVENTLIST which contains the eventflag to identify the types of artifact epoch has. As with manual visual inspection, marking/unmarking epochs will remove the identification of the type of artifact. It will just leave either the epoch has artifact(s) or not. So, there will be no (true) info about the kind(s) of artifact.
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 8,'setname',[filename, '-Epoched-ARTautomatic'],'savenew',[filepath, filesep, 'Marked Epochs - ', filter_type, filesep, filename, '-Epoched-ARTautomatic.set'],'gui','off');
            EEG.data = double(EEG.data);    % pop_newset with save, saves in single format as double precision is not required in saving.
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 9,'setname',[filename, '-Epoched-ARTmanual'],'savenew',[filepath, filesep, 'Marked Epochs - ', filter_type, filesep, filename, '-Epoched-ARTmanual.set'],'gui','off');
            EEG.data = double(EEG.data);    % pop_newset with save, saves in single format as double precision is not required in saving.
            EEG = pop_exporteegeventlist( EEG , 'Filename', [filepath, filesep, 'Marked Epochs - ', filter_type, filesep, filename, '-eventList-ARTautomatic.txt'] );
            
            [EEG, tprej, acce, rej, histoflags ] = pop_summary_AR_eeg_detection(EEG, [filepath, filesep, 'Marked Epochs - ', filter_type, filesep, filename, '-stats-ARTautomatic.txt']);   %prints table to file
            
            stats_ARTautomatic(length(stats_ARTautomatic)+1).subject = subject_id(i);
            % stats_ARTautomatic(length(stats_ARTautomatic)).session = str2double(filename);
            stats_ARTautomatic(length(stats_ARTautomatic)).session = filename;
            stats_ARTautomatic(length(stats_ARTautomatic)).accepted = acce;
            stats_ARTautomatic(length(stats_ARTautomatic)).rejected = rej;
            stats_ARTautomatic(length(stats_ARTautomatic)).percent_rejected = tprej;
            stats_ARTautomatic(length(stats_ARTautomatic)).F2_max_voltage_threshold = histoflags(:, 2);
            stats_ARTautomatic(length(stats_ARTautomatic)).F3_pp_threshold = histoflags(:, 3);
            stats_ARTautomatic(length(stats_ARTautomatic)).F4_step_like = histoflags(:, 4);
            stats_ARTautomatic(length(stats_ARTautomatic)).F5_sample_difference = histoflags(:, 5);
            stats_ARTautomatic(length(stats_ARTautomatic)).F6_flatline = histoflags(:, 6);
            
            %% mark/unmark epochs for rejection manually
            % 1) in EEGLAB gui click file > load existing dataset > session#-Epoched-ARTmanual.set
            % 2) run mark_epochs_gui.m code
            % 3) in mark_epochs_gui check all 6 boxes, 'Modified' & 'OR' then click 'plot for inspection' 
            % 4) mark/unmark epochs then when finished click 'update marks' > 'ok' > 'save (overwrites)' > 'yes' > close mark_epochs_gui
            % 5) repeat steps 1-4 with next -Epoched-ARTmanual file
        end
        
        
        %% ICA
        % mark epochs for rejection using ERPLAB, and verify manually if the
        % automatic detection of nosiy epochs worked fine. Mark/unmark epochs
        % manually, before running ICA.
        if (run_ICA)
            
            %% Get marked epochs for rejection
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Getting Epoch Marks for Rejection']);
            
            % check if file -Epoched-ARTmanual is loaded. Most likely it will
            % not be loaded, as after running epoch rejection, manual
            % inspection is done. Then, ICA is run. So, in almost all cases, it
            % will be loaded here.
            set_idx = [];
            if (~isempty(ALLEEG))
                tmp = {ALLEEG.setname};
                set_idx = find(strcmp(tmp, [filename, '-Epoched-ARTmanual']));
            end
            if (isempty(set_idx))
                load_filepath = [filepath, filesep, 'Marked Epochs - ', filter_type,];
                EEG = pop_loadset('filename',[filename, '-Epoched-ARTmanual.set'],'filepath',load_filepath);
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                EEG = eeg_checkset( EEG );
            else
                EEG = pop_newset(ALLEEG, EEG, CURRENTSET, 'retrieve', set_idx);
                EEG = eeg_checkset( EEG );
            end
            
            reject_epochs = EEG.reject.rejmanual;
            %         reject_epochs = tmp_reject;
            
            
            %% downsample continuous HighPass filtered file
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Downsampling']);
            
            % check if file -HPfiltered is loaded. Most likely it will
            % not be loaded, as after running epoch rejection, manual
            % inspection is done. Then, ICA is run. So, in almost all cases, it
            % will be loaded here.
            set_idx = [];
            if (~isempty(ALLEEG))
                tmp = {ALLEEG.setname};
                set_idx = find(strcmp(tmp, [filename, '-HPfiltered']));
            end
            if (isempty(set_idx))
                load_filepath = [filepath, filesep, 'Marked Epochs - ', filter_type,];
                EEG = pop_loadset('filename',[filename, '-HPfiltered.set'],'filepath',load_filepath);
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                EEG = eeg_checkset( EEG );
            else
                EEG = pop_newset(ALLEEG, EEG, CURRENTSET, 'retrieve', set_idx);
                EEG = eeg_checkset( EEG );
            end
            
            % downsample continuous file, as downsampling epoched data is not
            % recommended due to anti-aliasing filtering
            EEG = pop_overwritevent(EEG, 'code');   % added because after downsampling, ERPLAB doesn't work correctly with EEG.event.type having codelabels. The pop_editeventlist works OK then, and when finished returns EEG.event.type fileld with codelabels.
            
            EEG = pop_resample(EEG, ICA_downsample_rate);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 10,'setname',[filename, '-downsampled'],'gui','off');
            
            
            %% extract epochs
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Extracting Epochs']);
            
            % extract epochs from downsampled data
            EEG  = pop_editeventlist(EEG, ...
                'AlphanumericCleaning', 'on', ...
                'BoundaryNumeric', { -99}, ...
                'BoundaryString', { 'boundary' }, ...
                'List', [data_folder, 'data_info', filesep, 'elist-equations-', cell2mat(eeg_type), '.txt'], ...
                'SendEL2', 'EEG', ...
                'UpdateEEG', 'codelabel', ...
                'Warning', 'off' );
            
            if (strcmp (eeg_type, 'rest'))
                EEG = pop_epochbin(EEG, [0.0 1000*(EEG.etc.epoch_length-1/(2*EEG.srate))], 'none');
            elseif (strcmp (eeg_type, 'seps'))
                EEG = pop_epochbin(EEG, [-100.0 150.0], 'pre');     
            end
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 11,'setname',[filename, '-Epoched-Downsampled'],'gui','off');
            
            
            %% reject epochs
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Rejecting Epochs']);
            
            % reject marked epochs
            EEG = pop_rejepoch(EEG, reject_epochs, 0);
            % EEG = pop_syncroartifacts(EEG, 'Direction', 'bidirectional'); %
            % not required, as EEG dataset is not stored. Also bidirectional
            % only marks any unmarked epochs, and doesn't unmark any. So, not
            % correct in this case.
            
            
            %% ICA
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'ICA', '-started:', char(datetime)]);
            
            if (ICA_with_PCA)
                referenceChannels = EEG.etc.noiseDetection.reference.referenceChannels;
                interpolatedChannels = EEG.etc.noiseDetection.interpolatedChannelNumbers;
                channelsLeft = setdiff(referenceChannels, interpolatedChannels);
                pcaDim = length(channelsLeft) - 1;      % -1 because avg referenced
                
                [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 12,'setname',[filename, '-Ready-for-ICA'],'savenew',[ICA_temp_folder, filesep, filename, '-Ready-for-ICA'],'gui','off');
                EEG.data = double(EEG.data);    % pop_newset with save, saves in single format as double precision is not required in saving.
                EEG = eeg_checkset( EEG );
                
                runamica15(EEG, 'outdir', [ICA_temp_folder, filesep, 'amicaout'], ...
                    'pcakeep', pcaDim, 'max_threads', 2);      % aalborg crusher has processor with 2 threads
            else
                interpolatedChannels = EEG.etc.noiseDetection.interpolatedChannelNumbers;
                stillNoisyChannels = EEG.etc.noiseDetection.stillNoisyChannelNumbers;
                remove_channels = [interpolatedChannels, stillNoisyChannels];
                
                if (isempty(remove_channels))   % account for average referencing by removing 1 channel
                    remove_channels = {'POZ'};
                end
                
                EEG = pop_select(EEG, 'nochannel', remove_channels);
                [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 12,'setname',[filename, '-Ready-for-ICA'],'savenew',[ICA_temp_folder, filesep, filename, '-Ready-for-ICA'],'gui','off');
                EEG.data = double(EEG.data);    % pop_newset with save, saves in single format as double precision is not required in saving.
                EEG = eeg_checkset( EEG );
                
                runamica15(EEG, 'outdir', [ICA_temp_folder, filesep, 'amicaout'], ...
                    'pcakeep', EEG.nbchan, 'max_threads', 4);      % NZ computer with 4 physical cores & no hyperthreadng
            end
            
            % move file from ICA temp folder to the subject's folder, and
            % remove files from ICA temp folder
            previousState = recycle('off');     % permanently delete the files, instead of putting in the recycle bin
            movefile ([ICA_temp_folder, filesep, 'amicaout'], [filepath, filesep, 'amicaout - ', filter_type]);
            delete ([ICA_temp_folder, filesep, filename, '-Ready-for-ICA.set']);
            delete ([ICA_temp_folder, filesep, filename, '-Ready-for-ICA.fdt']);
            recycle(previousState);     % reset the settings to the previous one/default (can be set in the Preferences->General)
        end
        
        
        %% Apply IC matrix, calculate dipfit, SASICA (with ADJUST and FASTER)
        % mark ICs for rejection using mark_ICs_gui after applying weights and
        % calculating the required material here which helps in marking ICs
        if (apply_IC_weights_rejection)
            
            %% Get marked epochs for rejection
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Getting Epoch Marks for Rejection']);
            
            % check if file -Epoched-ARTmanual is loaded. Most likely it will
            % not be loaded, as after running epoch rejection, manual
            % inspection is done. Then, ICA is run. So, in almost all cases, it
            % will be loaded here.
            set_idx = [];
            if (~isempty(ALLEEG))
                tmp = {ALLEEG.setname};
                set_idx = find(strcmp(tmp, [filename, '-Epoched-ARTmanual']));
            end
            if (isempty(set_idx))
                load_filepath = [filepath, filesep, 'Marked Epochs - ', filter_type,];
                EEG = pop_loadset('filename',[filename, '-Epoched-ARTmanual.set'],'filepath',load_filepath);
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                EEG = eeg_checkset( EEG );
            else
                EEG = pop_newset(ALLEEG, EEG, CURRENTSET, 'retrieve', set_idx);
                EEG = eeg_checkset( EEG );
            end
            
            reject_epochs = EEG.reject.rejmanual;
            %         reject_epochs = tmp_reject;
            
            
            %% filter data
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Bandpass Filtering']);
            
            % check if file "with event list" is loaded, e.g. if PREP
            % pipeline was run earlier, and event list created,
            % and no longer required to repeat and is set to false,
            % then PREP file "with event list" will not be created and not
            % present in current datasets.
            set_idx = [];
            if (~isempty(ALLEEG))
                tmp = {ALLEEG.setname};
                set_idx = find(strcmp(tmp, [filename, '-eventListed']));
            end
            if (isempty(set_idx))
                load_filepath = [filepath, filesep, 'eventListed'];
                EEG = pop_loadset('filename',[filename, '-eventListed.set'],'filepath',load_filepath);
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                EEG = eeg_checkset( EEG );
            else
                EEG = pop_newset(ALLEEG, EEG, CURRENTSET, 'retrieve', set_idx);
                EEG = eeg_checkset( EEG );
            end
            
            if (strcmp(filter_type, 'FIR'))
                % the order of the filter was calculated using pop_firwsord, or
                % EEGLAB->Tools->Filter the data->Windowed Sinc FIR Filter
                EEG = pop_firws(EEG, 'fcutoff', [filter_bandpass.HP_cutoff, filter_bandpass.LP_cutoff], 'ftype', 'bandpass', 'wtype', 'kaiser', 'warg', 5.65326, 'forder', filter_order, 'minphase', 0);
            elseif (strcmp (filter_type, 'IIR'))
                [b, a] = butter (filter_order, [filter_bandpass.HP_cutoff, filter_bandpass.LP_cutoff]/(EEG.srate/2), 'bandpass');
                EEG.data = filtfilt (b, a, EEG.data')';
            end
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 13,'setname',[filename, '-filtered'],'gui','off');
            EEG = eeg_checkset( EEG );
            
            
            %% extract epochs
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Extracting Epochs']);
            
            
            
            if (strcmp (eeg_type, 'rest'))
                EEG = pop_epochbin(EEG, [0.0 1000*(EEG.etc.epoch_length-1/(2*EEG.srate))], 'none');
            elseif (strcmp (eeg_type, 'seps'))
                EEG = pop_epochbin(EEG, [-100.0 150.0], 'pre');     
            end
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 14,'setname',[filename, '-Epoched'],'gui','off');
            
            
            %% reject epochs
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Rejecting Epochs']);
            
            % reject marked epochs
            EEG = pop_rejepoch(EEG, reject_epochs, 0);
            % EEG = pop_syncroartifacts(EEG, 'Direction', 'bidirectional'); %
            % not required, as EEG dataset is not stored. Also bidirectional
            % only marks any unmarked epochs, and doesn't unmark any. So, not
            % correct in this case.
            
            
            %% apply ICA weights
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Applying ICA weights']);
            
            EEG.etc.chanlocs = EEG.chanlocs;    % save channel locations to use in interpolation after rejection of ICs
            
            filepath_amica = [filepath, filesep, 'amicaout - ', filter_type,];
            
            if (ICA_with_PCA)       % need to fix this ICA with PCA by verifying how many components are returned when ICA is done with PCA with avg reference, and how to apply ICA weights to EEG data.
                             referenceChannels = EEG.etc.noiseDetection.reference.referenceChannels;
                             interpolatedChannels = EEG.etc.noiseDetection.interpolatedChannelNumbers;
                             channelsLeft = setdiff(referenceChannels, interpolatedChannels);
                             pcaDim = length(channelsLeft) - 1;      % -1 because avg referenced
                
                             [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 12,'setname',[filename, '-Epoched-Rejected'],'savenew',[ICA_temp_folder, filesep, filename, '-Epoched-Rejected'],'gui','off');
                             EEG = eeg_checkset( EEG );
                
                             runamica15(EEG, 'outdir', [ICA_temp_folder, filesep, 'amicaout'], ...
                                 'pcakeep', pcaDim, 'max_threads', 2);      % aalborg crusher has processor with 2 threads
            else
                interpolatedChannels = EEG.etc.noiseDetection.interpolatedChannelNumbers;
                stillNoisyChannels = EEG.etc.noiseDetection.stillNoisyChannelNumbers;
                remove_channels = [interpolatedChannels, stillNoisyChannels];
                
                if (isempty(remove_channels))   % account for average referencing by removing 1 channel
                    remove_channels = {'POZ'};
                end
                
                EEG = pop_select(EEG, 'nochannel', remove_channels);
            end
            
            save_filepath = [filepath, filesep, 'ICA - ', filter_type];
            if (~exist(save_filepath, 'dir'))
                mkdir (save_filepath);
            end
            
            %% apply ICA matrix
            EEG = eeg_loadamica(EEG, filepath_amica, 1);
            
                       
            %%% run dipfit with single, bilateral and fitTwoDipoles dipole fitting
            %% to get the code for EEG = pop_dipfit_settings(...) below follow these steps...
            %% immediately following selection of channel locations above, apply ICA weights on Matlab command line type EEG.icasphere = eye(EEG.nbchan); > press Enter > EEG.icaweights = eye(EEG.nbchan); > press Enter > eeglab redraw > Enter
            %% in EEGLAB GUI ICA weights should now be labeled Yes then select Tools > Locate dipoles using Dipfit > select Boundary Element Model (MNI) > select channel to omit (if any) > Manual Co-Reg > resize X,Y,Z to 1.5 > select green & brown Labels On > select brown Electrodes (21 elec, 10-20) > Warp montage > Ok > Ok > Ok 
            %% on Matlab command line type eegh
            %EEG = pop_dipfit_settings( EEG, 'hdmfile','E:\\Folders\\PhD\\ISAD\\ISAD_Curry_data\\ISAD_EEGs\\Preprocessing\\eeg_preprocessing\\eeglab14_1_1b\\eeglab14_1_1b\\plugins\\dipfit2.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','E:\\Folders\\PhD\\ISAD\\ISAD_Curry_data\\ISAD_EEGs\\Preprocessing\\eeg_preprocessing\\eeglab14_1_1b\\eeglab14_1_1b\\plugins\\dipfit2.3\\standard_BEM\\standard_mri.mat','chanfile','E:\\Folders\\PhD\\ISAD\\ISAD_Curry_data\\ISAD_EEGs\\Preprocessing\\eeg_preprocessing\\eeglab14_1_1b\\eeglab14_1_1b\\plugins\\dipfit2.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[0.96125 -17.6702 -1.1575 -3.2029e-07 3.4264e-07 -1.5708 1 1 1] ,'chansel',[1:60] );
            EEG = pop_dipfit_settings( EEG, 'hdmfile','E:\\Folders\\PhD\\ISAD\\ISAD_Curry_data\\ISAD_EEGs\\Preprocessing\\eeg_preprocessing\\eeglab14_1_1b\\eeglab14_1_1b\\plugins\\dipfit2.3\\standard_BEM\\standard_vol.mat','coordformat','MNI','mrifile','E:\\Folders\\PhD\\ISAD\\ISAD_Curry_data\\ISAD_EEGs\\Preprocessing\\eeg_preprocessing\\eeglab14_1_1b\\eeglab14_1_1b\\plugins\\dipfit2.3\\standard_BEM\\standard_mri.mat','chanfile','E:\\Folders\\PhD\\ISAD\\ISAD_Curry_data\\ISAD_EEGs\\Preprocessing\\eeg_preprocessing\\eeglab14_1_1b\\eeglab14_1_1b\\plugins\\dipfit2.3\\standard_BEM\\elec\\standard_1005.elc','coord_transform',[0.96125 -17.6702 -1.1575 -3.2029e-07 3.4264e-07 -1.5708 1 1 1] ,'chansel',[1:EEG.nbchan] );


            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
 
            % single dipole
            EEG = pop_multifit(EEG, [1:EEG.nbchan], 'threshold', 100, 'plotopt', {'normlen' 'on'});
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
            EEG.etc.dipfit.single_dipole = EEG.dipfit;
            
            % bilateral dipoles
            EEG = pop_multifit(EEG, [1:EEG.nbchan], 'threshold', 100, 'dipoles', 2, 'plotopt', {'normlen' 'on'});
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
            EEG.etc.dipfit.bilateral_dipole = EEG.dipfit;
            
            % fitTwoDipoles
            % pop_dipplot( EEG,[1:57] ,'mri','/Applications/MATLAB_R2015b.app/toolbox/eeglab14_1_0b/plugins/dipfit2.3/standard_BESA/avg152t1.mat','normlen','on');
            try         %fitTwoDipoles fails when I think it does not converge, and so rv variable is not generated and code crashes
                EEG = fitTwoDipoles(EEG, 'LRR', 35);
                
                EEG.etc.dipfit.fitTwo_dipole_success = true;
                EEG.etc.dipfit.fitTwo_dipole = EEG.dipfit;
            catch matlabException
                EEG.etc.dipfit.fitTwo_dipole_success = false;
                
                errordlg({['Error fitTwoDipoles in file - ', filename],...
                    ['Identifier: ', matlabException.identifier], ...
                    matlabException.message}, ...
                    'Error', 'non-modal');
            end
            
            % set EEG dipfit struct back to single dipole dipfit
            EEG.dipfit = EEG.etc.dipfit.single_dipole;
            
            %% 'percent variance accounted for' (pvaf) and  power spectral density (PSD)
            
            for idx = 1:length(EEG.icachansind)
                EEG.etc.ica_pvaf(idx, 1) = eeg_pvaf (EEG, idx, 'plot', 'off');
                
                EEG.etc.ica_psd_spectra(idx, :) = spectopo (EEG.icaact(idx,:), EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,idx), 'nfft', EEG.srate, 'plot', 'off');        %freqfac = 4 for 2048/512 gives 1Hz frequency resolution of fft
            end
            [~, EEG.etc.ica_psd_freq] = spectopo (EEG.icaact(idx,:), EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,idx), 'nfft', EEG.srate, 'plot', 'off');
            
            %% SASICA
            % find the available channels which can be given as blink channels
            % for FASTER. note: assuming at least 1 of the 'FP1', 'FPz', 'FP2',
            % 'AF7', 'AF3', 'AF4', 'AF8' exists. If none of 'FP1', 'FPz', 'FP2'
            % exists then look for 'AF7', 'AF3', 'AF4', 'AF8'. Otherwise don't
            % do FASTER algorithm.
            enable_FASTER = 1;
            enableADJUST = 1;
            blink_channels_1 = {'FP1', 'FPz', 'FP2'};
            blink_channels_2 = {'AF7', 'AF3', 'AF4', 'AF8'};
            chan_labels = {EEG.chanlocs.labels};
            blink_channels = ismember (chan_labels, blink_channels_1);
            if (sum(blink_channels) < 1)
                blink_channels = ismember (chan_labels, blink_channels_2);
                
                if (sum(blink_channels) < 1)
                    enable_FASTER = 0;
                    enableADJUST = 0;
                end
            end
            blink_channels = chan_labels (blink_channels);
            
            [EEG] = eeg_SASICA (EEG,'MARA_enable',0, 'FASTER_enable',enable_FASTER,'FASTER_blinkchanname',blink_channels, ...
                'ADJUST_enable',enableADJUST,'chancorr_enable',0,'chancorr_channames','No channel', ...
                'chancorr_corthresh','auto 4','EOGcorr_enable',0,'EOGcorr_Heogchannames','No channel', ...
                'EOGcorr_corthreshH','auto 4','EOGcorr_Veogchannames','No channel','EOGcorr_corthreshV','auto 4', ...
                'resvar_enable',1,'resvar_thresh',15,'SNR_enable',0,'SNR_snrcut',1,'SNR_snrBL',[-Inf 0], ...
                'SNR_snrPOI',[0 Inf] ,'trialfoc_enable',1,'trialfoc_focaltrialout','auto', ...
                'focalcomp_enable',1,'focalcomp_focalICAout','auto','autocorr_enable',1,'autocorr_autocorrint',20, ...
                'autocorr_dropautocorr','auto','opts_noplot',1,'opts_nocompute',0,'opts_FontSize',14);
            
            %% save
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 15,'setname',[filename, '-ICA'],'savenew',[filepath, filesep, 'ICA - ', filter_type, filesep, filename, '-ICA.set'],'gui','off');
            EEG.data = double(EEG.data);    % pop_newset with save, saves in single format as double precision is not required in saving.
        end
        
        
        %% reject ICA components
        % 1) in EEGLAB gui click 'file' > 'load existing dataset' & select desired -ICA.set file
        % 2) run mark_ICs_gui.m code
        % 3) label ICs as 'Brain' or 'Other' (only necessary for first ~30 ICs) then right click on '?' to mark remaining ICs as 'Others') 
        % 4) click 'save (overwrites)' > 'yes'
        % 5) visualise plot of blue uncleaned &  red cleaned (i.e., 'Other' ICs removed)
        % 6) can double-check rejected ICs in EEGLAB gui by clicking 'Tools' > 'Reject data using ICA' > 'Reject components by map'(green kept, red removed)
        % 7) repeat steps 1-5 for next -ICA.set file
        % 8) set reject_ICA_components to 'true'
  
        if (reject_ICA_components)
            
            %% Get marked epochs for updating flags and EEG.reject struct
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Getting Epoch Marks for Updating Flags']);
            
            % check if file -Epoched-ARTmanual is loaded. Most likely it will
            % not be loaded, as after running epoch rejection, manual
            % inspection is done. Then, ICA is run. So, in almost all cases, it
            % will be loaded here.
            set_idx = [];
            if (~isempty(ALLEEG))
                tmp = {ALLEEG.setname};
                set_idx = find(strcmp(tmp, [filename, '-Epoched-ARTmanual']));
            end
            if (isempty(set_idx))
                load_filepath = [filepath, filesep, 'Marked Epochs - ', filter_type,];
                EEG = pop_loadset('filename',[filename, '-Epoched-ARTmanual.set'],'filepath',load_filepath);
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                EEG = eeg_checkset( EEG );
            else
                EEG = pop_newset(ALLEEG, EEG, CURRENTSET, 'retrieve', set_idx);
                EEG = eeg_checkset( EEG );
            end
            
            reject_epochs.rejmanual = EEG.reject.rejmanual;
            reject_epochs.rejmanualE = EEG.reject.rejmanualE;
            
            eventlist_flags = {EEG.EVENTLIST.eventinfo.flag};
            
            
            %% Get marked ICs
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Getting Marked ICs']);
            
            % check if file -ICA is loaded. Most likely it will
            % not be loaded, as after running ICA manual
            % inspection of ICs is done. So, in almost all cases, it
            % will be loaded here.
            set_idx = [];
            if (~isempty(ALLEEG))
                tmp = {ALLEEG.setname};
                set_idx = find(strcmp(tmp, [filename, '-ICA']));
            end
            if (isempty(set_idx))
                load_filepath = [filepath, filesep, 'ICA - ', filter_type,];
                EEG = pop_loadset('filename',[filename, '-ICA.set'],'filepath',load_filepath);
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                EEG = eeg_checkset( EEG );
            else
                EEG = pop_newset(ALLEEG, EEG, CURRENTSET, 'retrieve', set_idx);
                EEG = eeg_checkset( EEG );
            end
            
            % marked/unmarked IC components
            reject_IC_components = EEG.reject.gcompreject;
            EEG_chanlocs = EEG.etc.chanlocs;            % channel locations saved when ICA matrix was applied to reject the ICs, retrieved here for interpolation after applying ICA matrix
            
            
            %% filter data
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Bandpass Filtering']);
            
            % check if file "with event list" is loaded, e.g. if PREP
            % pipeline was run earlier, and event list created,
            % and no longer required to repeat and is set to false,
            % then PREP file "with event list" will not be created and not
            % present in current datasets.
            set_idx = [];
            if (~isempty(ALLEEG))
                tmp = {ALLEEG.setname};
                set_idx = find(strcmp(tmp, [filename, '-eventListed']));
            end
            if (isempty(set_idx))
                load_filepath = [filepath, filesep, 'eventListed'];
                EEG = pop_loadset('filename',[filename, '-eventListed.set'],'filepath',load_filepath);
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                EEG = eeg_checkset( EEG );
            else
                EEG = pop_newset(ALLEEG, EEG, CURRENTSET, 'retrieve', set_idx);
                EEG = eeg_checkset( EEG );
            end
            
%             if (strcmp(filter_type, 'FIR'))
%                 % the order of the filter was calculated using pop_firwsord, or
%                 % EEGLAB->Tools->Filter the data->Windowed Sinc FIR Filter
%                 EEG = pop_firws(EEG, 'fcutoff', [filter_bandpass.HP_cutoff, filter_bandpass.LP_cutoff], 'ftype', 'bandpass', 'wtype', 'kaiser', 'warg', 5.65326, 'forder', filter_order, 'minphase', 0);
%             elseif (strcmp (filter_type, 'IIR'))
                [b, a] = butter (filter_order, [filter_bandpass.HP_cutoff, filter_bandpass.LP_cutoff]/(EEG.srate/2), 'bandpass');
                EEG.data = filtfilt (b, a, EEG.data')';
%             end
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 16,'setname',[filename, '-filtered'],'gui','off');
            EEG = eeg_checkset( EEG );
            
            
            %% extract epochs
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Extracting Epochs']);
            
            if (strcmp (eeg_type, 'rest'))
                EEG = pop_epochbin(EEG, [0.0 1000*(EEG.etc.epoch_length-1/(2*EEG.srate))], 'none');
            elseif (strcmp (eeg_type, 'seps'))
                EEG = pop_epochbin(EEG, [-100.0 150.0], 'pre');     
            end
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 17,'setname',[filename, '-Epoched'],'gui','off');
            
            %% cleaned data without ICA
            % data without artifact epochs. although the epochs are not
            % rejected but the information of bad epochs is retained, so it
            % can be used in ERPLAB if required.
            
            % update EVENTLIST, EEG.epoch flag struct and
            % EEG.reject.rejmanual
            for flg = 1:EEG.trials
                EEG.EVENTLIST.eventinfo(flg).flag = eventlist_flags{flg};
                EEG.epoch(flg).eventflag = eventlist_flags{flg};
            end
            EEG.reject.rejmanual = reject_epochs.rejmanual;
            EEG.reject.rejmanualE = reject_epochs.rejmanualE;
            EEG = eeg_checkset( EEG );
            
%           save_filepath = [filepath, filesep, 'Epochs-cleaned - ', filter_type, '(', num2str(filter_bandpass.HP_cutoff), '-', num2str(filter_bandpass.LP_cutoff), ')'];
            save_filepath = [filepath, filesep, 'Epochs-cleaned - ', 'IIR', '(', num2str(filter_bandpass.HP_cutoff), '-', num2str(filter_bandpass.LP_cutoff), ')'];

            if (~exist(save_filepath, 'dir'))
                mkdir (save_filepath);
            end
            
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 17,'setname',[filename, '-Epochs-cleaned'],'savenew',[save_filepath, filesep, filename, '-Epochs-cleaned.set'],'gui','off');
            EEG.data = double(EEG.data);    % pop_newset with save, saves in single format as double precision is not required in saving.
            
            %% apply ICA weights
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Applying ICA weights']);
            
            %         EEG.etc.chanlocs = EEG.chanlocs;    % save channel locations to use in interpolation after rejection of ICs
            
            filepath_amica = [filepath, filesep, 'amicaout - ', filter_type];
            
            if (ICA_with_PCA)       % need to fix this ICA with PCA by verifying how many components are returned when ICA is done with PCA with avg reference, and how to apply ICA weights to EEG data.
                %             referenceChannels = EEG.etc.noiseDetection.reference.referenceChannels;
                %             interpolatedChannels = EEG.etc.noiseDetection.interpolatedChannelNumbers;
                %             channelsLeft = setdiff(referenceChannels, interpolatedChannels);
                %             pcaDim = length(channelsLeft) - 1;      % -1 because avg referenced
                %
                %             [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 12,'setname',[filename, '-Epoched-Rejected'],'savenew',[ICA_temp_folder, filesep, filename, '-Epoched-Rejected'],'gui','off');
                %             EEG = eeg_checkset( EEG );
                %
                %             runamica15(EEG, 'outdir', [ICA_temp_folder, filesep, 'amicaout'], ...
                %                 'pcakeep', pcaDim, 'max_threads', 2);      % aalborg crusher has processor with 2 threads
            else
                interpolatedChannels = EEG.etc.noiseDetection.interpolatedChannelNumbers;
                stillNoisyChannels = EEG.etc.noiseDetection.stillNoisyChannelNumbers;
                remove_channels = [interpolatedChannels, stillNoisyChannels];
                
                if (isempty(remove_channels))   % account for average referencing by removing 1 channel
                    remove_channels = {'POZ'};
                end
                
                EEG = pop_select(EEG, 'nochannel', remove_channels);
            end
            
            % the directory should exist, as it was created when ICA matrix was
            % applied
            %         load_filepath = [filepath, filesep, 'ICA - ', filter_type];
            %         if (~exist(load_filepath, 'dir'))
            %             mkdir (load_filepath);
            %         end
            
            %% apply ICA matrix
            EEG = eeg_loadamica(EEG, filepath_amica, 1);
            
            
            %% reject components
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Rejecting ICs']);
            
            % reject components
            EEG = pop_subcomp(EEG, find(reject_IC_components == 1), 0);
            
            % interpolation of channels
            EEG = pop_interp(EEG, eeg_mergelocs(EEG_chanlocs), 'spherical');
            
            % do baseline correction again for seps: this is done when plotting
            % SEPs using ERPLAB, and also in the code in the mark_SEPs_gui.
            
            
            %% update epoch flags
            % this helps in keeping all epochs with any marks, so ERPLAB can be
            % used with more power, e.g. plotting avg ERPs having all ERPs,
            % only clean ERPs or only artifact/noisy ERPs
            % done this at the end, as pop_select resets the EEG.epoch struct
            % and , pop_subcomp, pop_interp resets the EEG.reject(.rejmanual)
            % I think be careful to use EEGLAB functions on this dataset.
            % ERPLAB will likely work fine. If not, then may be make a custom
            % function to synchronize these after every EEGLAB function which
            % cause issues.
            waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, ...
                [num2str(i), '/', num2str(length(subject_id)), '-', 'Updating Epoch Flags']);
            
            % update EVENTLIST, EEG.epoch flag struct and
            % EEG.reject.rejmanual
            for flg = 1:EEG.trials
                EEG.EVENTLIST.eventinfo(flg).flag = eventlist_flags{flg};
                EEG.epoch(flg).eventflag = eventlist_flags{flg};
            end
            EEG.reject.rejmanual = reject_epochs.rejmanual;
            EEG.reject.rejmanualE = reject_epochs.rejmanualE;
            EEG = eeg_checkset( EEG );
            
            %% save
%           save_filepath = [filepath, filesep, 'ICA-cleaned - ', filter_type, '(', num2str(filter_bandpass.HP_cutoff), '-', num2str(filter_bandpass.LP_cutoff), ')'];
            save_filepath = [filepath, filesep, 'ICA-cleaned - ', 'IIR', '(', num2str(filter_bandpass.HP_cutoff), '-', num2str(filter_bandpass.LP_cutoff), ')'];

            if (~exist(save_filepath, 'dir'))
                mkdir (save_filepath);
            end
            
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 18,'setname',[filename, '-ICA-cleaned'],'savenew',[save_filepath, filesep, filename, '-ICA-cleaned.set'],'gui','off');
            EEG.data = double(EEG.data);    % pop_newset with save, saves in single format as double precision is not required in saving.
        end
        
    catch matlabException
        missing_files(length(missing_files)+1).subject = num2str(data_info(i).subject);
        missing_files(length(missing_files)).session = data_info(i).session;
        missing_files(length(missing_files)).error = matlabException.message;
    end
    
    %% reset the variables for next iteration/file
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
end


if (checking_false_triggers)
    save ([data_folder, 'data_info', filesep, 'epoch_info_automatic_to_be_verified.mat'], 'epoch_info');
end

if (mark_epochs_rejection)
    stats_file = [data_folder, 'data_info', filesep, 'epoch_artifacts_stats_', filter_type, '.mat'];
    % append: serves 2 purposes.
    % 1: when not doing all subjects at once, it concatenates the subgroup
    % of subjects used to previous subset of subjects.
    % 2: running the same subgroup of subjects for checking different
    % parameters for artifact detection, they are concatenated for easier
    % comparison to see the performance of different parameters e.g. by
    % seeing how number of epochs marked as artifacts is changed and which
    % epochs are marked/unmarked when eegploted.
    % Delete the file after finalizing the parameters, and re-run the code
    % to only save one instance of each subject.
    if (exist (stats_file, 'file'))
        save (stats_file, 'stats_ARTautomatic', '-append');
    else
        save (stats_file, 'stats_ARTautomatic');
    end
    
    
    rejection_parameters_file = [data_folder, 'data_info', filesep, 'rejection_parameters_', filter_type, '.mat'];
    
    % appends to existing file if already exists for comparion purposes
    if (exist (rejection_parameters_file, 'file'))
        save (rejection_parameters_file, 'epoch_rej_parameters', '-append');
    else
        save (rejection_parameters_file, 'epoch_rej_parameters');
    end
    
end

save ([data_folder, 'data_info', filesep, 'missing_files_ICA.mat'], 'missing_files');

%close EEGLAB gui
h = findobj('tag','EEGLAB');
close (h);

close(w);
% clc

% clear

rmpath ('E:\Folders\PhD\ISAD\ISAD_Curry_data\ISAD_EEGs\Preprocessing\eeg_preprocessing\eeglab14_1_1b\eeglab14_1_1b')

end