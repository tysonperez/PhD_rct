% montage2eeglab
% INPUT: montage - gtec Montage used for measuring the positions
% positionMatrix - measured postions in format (channels x 3 for x,y,z)
% filename - filename to store the converted positions

function montage2EEGlab(montage,filename)
    % load information from montage-set
    % check montage to be of montage type
    if ~strcmpi(class(montage),'montage')
        errordlg('Invalid montage file loaded!gtec montage file created by g.MONcreator is expected!','Invalid file');
        return;
    end
    channels = montage(11);
    numChannels = size(channels, 2);
    nameChannels = montage(12);
    xPosition = montage(16);
    yPosition = montage(17);
    zPosition = montage(18);
    
    % sort channels, channelnames (ascending channels)
    [channels, index] = sort(channels);
    nameChannels = nameChannels(index);
    xPosition = xPosition(index);
    yPosition = yPosition(index);
    zPosition = zPosition(index);

    file = fopen(filename,'w');
    for i=1:numChannels
        fprintf('%d %.4f %.4f %.4f %s\n',i,xPosition(i),yPosition(i),zPosition(i),nameChannels{i});
        fprintf(file,'%d %.4f %.4f %.4f %s\n',i,xPosition(i),yPosition(i),zPosition(i),nameChannels{i});
    end
    fclose(file);
end