% gtec import function

function [EEGout, command] = pop_read_gtecposition(EEG, filename)

    command = 'pop_read_gtecposition(ALLEEG)';

    fprintf('g.tec location file import function');

    if nargin < 3
        % ask user for file to import
        % may we include *.hdf5 format
        [filenametmp, filepath] = uigetfile('*.xyz', 'Choose a g.tec location file -- pop_read_gtecposition');
        drawnow;
        if filenametmp == 0 return; end;
        filename = [filepath '/' filenametmp];
    end;

    %locations = load(filename);

    fprintf('g.tec location file selected');

    fprintf('pop_read_gtecposition: importing locations from gtec montage file...\n');
    
%% input of montage creator file
%     montage = load(filename);
%     montage = montage.Mon;
%     if ~isobject(montage)
%         return;
%     end;
% 
%     channels = montage(11);
%     numChannels = size(channels, 2);
%     nameChannels = montage(12);
%     
%     positionMatrix = zeros(numChannels,3);
%     
%     positionMatrix(:,1) = montage(16)';
%     positionMatrix(:,2) = montage(17)';
%     positionMatrix(:,3) = montage(18)';
% 
%     % sort channels, channelnames (ascending channels)
%     [channels, index] = sort(channels);
%     nameChannels = nameChannels(index);
% 
%     file = fopen('gtec_tmp.xyz','w');
%     for i=1:numChannels
%         fprintf('%d %.4f %.4f %.4f %s\n',i,positionMatrix(i,1),positionMatrix(i,2),positionMatrix(i,3),nameChannels{i});
%         fprintf(file,'%d %.4f %.4f %.4f %s\n',i,positionMatrix(i,1),positionMatrix(i,2),positionMatrix(i,3),nameChannels{i});
%     end
%     fclose(file);
    
    EEGout = pop_chanedit(EEG, 'load', {filename 'filetype' 'xyz'}); 
    % save epochs in new datasets
    %[EEGout] = pop_editset(EEG, 'chanlocs', EEG.chanlocs);
    %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG_targets, 1, 'setname', [datasetname ' ' 'targetepochs']);
    %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG_nontargets, 1, 'setname', [datasetname ' ' 'nontargetepochs']);
        

    delete gtec_tmp.xyz;
return;
