%------------------------------------------------
%IMPORT CURRY (.cdt) FILES TO EEGLAB
% ------------------------------------------------
%%% folders should be set up Raw_EEGs > Suject1, Subject2,... > 1,2,3 > raw
%%% files (e.g. *.cdt)

save_filepath = ['Test'];

data_path = 'E:\Folders\PhD\ISAD\ISAD_Curry_data\ISAD_EEGs\Test';

missing_files = [];     % to store info of datasets not converted to eeglab format or missing files

save_filepath2 = [save_filepath, filesep, 'data_info '];
if (~exist(save_filepath2, 'dir'))
    mkdir (save_filepath2);
end

data_info = [];
di = 1;

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Get a list of all files and folders in this folder.
files = dir(data_path)
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir]
% Extract only those that are directories.
subFolders = files(dirFlags)                

% Print folder names to command window.
for folderi = 3:length(subFolders)
    fprintf('Sub folder #%d = %s\n', folderi, subFolders(folderi).name);
    
    files2 = dir([subFolders(folderi).folder, filesep, subFolders(folderi).name]);     
    
    dirFlags2 = [files2.isdir]
    
    % Extract only those that are directories.
    subFolders2 = files2(dirFlags2)                
    
    % loop through all folders 1,2,3,4,5
    for sessionFolderi = 3:length(subFolders2)
        files3 = dir([subFolders2(sessionFolderi).folder, filesep, subFolders2(sessionFolderi).name]);
        
        dirInfo = dir([files3(1).folder, filesep, '*.cdt']);
        
        filename = dirInfo.name;
        %             filename = convertCharsToStrings(filename);
        
        directory = [dirInfo.folder, filesep];
        %             directory = convertCharsToStrings(directory);
        %             directory = strcat(directory, "\", convertCharsToStrings(dirInfo(FN).name));
        %             directory = char(directory);
        
        % you need to fix this conditional statements.
        if (sessionFolderi == 3)
            session = 'session1';
        elseif (sessionFolderi == 4)
            session = 'session2';
        elseif (sessionFolderi == 5)
            session = 'session3';
        elseif (sessionFolderi == 6)
            session = 'session4';
        elseif (sessionFolderi == 7)
            session = 'session5';
        end
        
        % load file
        try
            EEG = loadcurry([directory, filename], 'CurryLocations', 'True', 'KeepTriggerChannel', 'False');
            EEG.setname = session;
            EEG = eeg_checkset( EEG );
            
            % make folder if it doesn't exist in the directory as
            % Data_EEGLAB/Subject 1/session name/Raw/
            % Data_EEGLAB/Subject 1/session2 name/Raw/
            % Data_EEGLAB/Subject 1/session3 name/Raw/ etc
            
            sub_number = files2(1).folder;
            sub_number = strsplit(sub_number, filesep);
            sub_number = sub_number{end};
            sub_number = strsplit(sub_number, 'H');
            sub_number = sub_number{end};
            
            save_filepath2 = [save_filepath, filesep, 'Subject ', sub_number, filesep, session, filesep, 'Raw'];
            if (~exist(save_filepath2, 'dir'))
                mkdir (save_filepath2);
            end
            
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname',[session, '-Raw'],'savenew',[save_filepath2, filesep, session, '-Raw.set'],'gui','off');
            
            % fill data_info
            data_info(di).session = session;
            data_info(di).subject = sub_number;
            data_info(di).setname = EEG.setname;
            data_info(di).comments = filename;
            data_info(di).EEG_length = EEG.xmax/60;
            di = di + 1;
            
            ALLEEG = pop_delset( ALLEEG, [1:length(ALLEEG)] );   %save memory
            
            fprintf('\n %s Import complete. \n \n', session)
            %eeglab redraw       % you don't need this to save time
            
            
        catch matlabException
            missing_files(length(missing_files)+1).subject = files2(1).folder;
            missing_files(length(missing_files)).session = session;
            missing_files(length(missing_files)).error = matlabException.message;
        end
    end
    
end

save ([save_filepath, filesep, 'data_info', filesep, 'missing_files_ICA.mat'], 'missing_files');
save ([save_filepath, filesep, 'data_info', filesep, 'data_share.mat'], 'data_info');

disp('Importing complete')
