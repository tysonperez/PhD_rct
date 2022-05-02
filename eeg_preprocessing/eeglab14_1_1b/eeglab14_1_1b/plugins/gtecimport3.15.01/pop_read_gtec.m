% gtec import function

function [EEG, ALLEEG, command] = pop_read_gtec(ALLEEG, filename)
EEG = [];
command = 'pop_read_gtec(ALLEEG)';

fprintf('g.tec file format import function');

if nargin < 2
	% ask user for file to import
    % may we include *.hdf5 format
	[filenametmp, filepath] = uigetfile('*.mat; *hdf5', 'Choose a g.tec file -- pop_read_gtec'); 
    drawnow;
	if filenametmp == 0 return; end;
	filename = [filepath  filenametmp];
end;

fprintf('g.tec file selected');

%  prompt user to insert sample rate, trigger window size and dataset name

% gtec supports 2 different file formats -> plain matlab array (*.mat) and
% extended format containing trigger infos, sampling rate,... (*.hdf)

% get which file was selected to be loaded

% read gtec format
% initialize EEG dataset
EEG = eeg_emptyset;
fprintf('pop_read_gtec: importing gtec file...\n');

[pathstr, name, ext] = fileparts(filename);
% set filename of dataset
EEG.filename = filename;

switch ext
    case '.mat'
        
        matDataArray=load(filename);
        
        % check input structure if timeseries object or matlab data array
        fieldNames = fieldnames(matDataArray);
        
        if(length(fieldNames) > 1)
            error('Invalid data file selected!');
        end
        
        field = getfield(matDataArray,fieldNames{1});
        % timeseries data or matlab array
        if isnumeric(field)
            % matlab array    
            sampleRate = 1 / (field(1,2) - field(1,1));
        elseif isa(field,'timeseries')
            % timeseries object
            sampleRate = 1 / (field.Time(2) - field.Time(1));
        else
            error('Invalid data structure')
        end
        
        
        prompt={'Check sampling rate of imported data:',...
            'Enter name of dataset:'...
            %'Enter epoch limits in relation to event:',...
            };
        name='Input for gtec import function';
        numlines=1;
        defaultanswer={sampleRate,'gtec_import'};
        %defaultanswer={'256','gtec_import','[-0.1 0.4]'};
        
        answer=inputdlg2(prompt,name,numlines,defaultanswer);
        
        for i=1:size(prompt,2)
            [val status] = str2num(answer{1});  % Use curly bracket for subscript
            if ~status
                % Handle empty value returned for unsuccessful conversion
                % ...
                %error('Unexpected input...');
            end
            
            % get input to variables
            switch i
                case 1
                    samplerate = str2double(answer{i});
                case 2
                    datasetname = answer{i};
                case 3
                    epochlimits = str2num(answer{i});
                otherwise
                    fprintf('Unhandled user input..');
                    return;
            end
        end

        % [EEG.data,events,header] = read_gtec(filename);
        % import data to EEG dataset
        EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',filename,'setname',datasetname,'srate',samplerate,'pnts',0,'xmin',0);
        EEG = eeg_checkset(EEG);
        
        load_classinfo = questdlg2('Does the data contain event information?','Event Information','Yes','No','Yes');

        if strcmpi(load_classinfo,'Yes')
            while 1
                prompt={'Enter event information channel number:'
                       % 'Enter name of event:',...
                };
                name='Input for event information';
                numlines=1;
                defaultanswer={EEG.nbchan};%,'TARGET'};
                %defaultanswer={'256','gtec_import','[-0.1 0.4]'};

                answer=inputdlg2(prompt,name,numlines,defaultanswer);
                if length(answer) ~= 1
                    errordlg2('Invalid event info configuration!','Event Information configuration');
                    break;
                end
                
                if isnumeric(answer{1}) && answer{1} < 1 && answer{1} > EEG.nbchan 
                    errordlg2('Invalid channel number entered!','Event Information channel');
                    break;
                end

                % import events from datachannel
                EEG = pop_chanevent(EEG,str2num(answer{1}),'edge','leading','edgelen',1,'delchan','off','delevent','off');
                EEG = eeg_checkset(EEG); 
                
                cont = questdlg2('Are there more channels with event information?','Event Information','Yes','No','Yes');
                if strcmpi(cont,'No')
                    break;
                end
            end
        else
        end
        
        % remove time channel
%         EEG.data([1],:) = [];
%         EEG.nbchan = EEG.nbchan - 1;
        
        % get event information
%         eventno = size(unique(EEG.data(EEG.nbchan-1,EEG.data(EEG.nbchan-1,:)~=0)),2)
%         targetevent = num2cell(EEG.data(EEG.nbchan-1,find(EEG.data(EEG.nbchan,:),1,'first')));
%         nontargetevent = num2cell(unique(EEG.data(EEG.nbchan-1,EEG.data(EEG.nbchan-1,:)~=9 & EEG.data(EEG.nbchan-1,:)~=0)));
%         
%         % create epochs
%         EEG_targets = pop_epoch(EEG,targetevent,epochlimits,'new_name',[datasetname ' ' 'targetepochs']);
%         EEG_nontargets = pop_epoch(EEG,nontargetevent,epochlimits,'new_name',[datasetname ' ' 'nontargetepochs']);
%         
%         % save epochs in new datasets
%         [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', datasetname);
%         [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG_targets, 1, 'setname', [datasetname ' ' 'targetepochs']);
%         [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG_nontargets, 1, 'setname', [datasetname ' ' 'nontargetepochs']);
    case '.hdf5'
        
        prompt={'Enter name of dataset:'...
            %'Enter epoch limits in relation to event:',...
            };
        name='Input for gtec import function';
        numlines=1;
        defaultanswer={'gtec_import'};
        %defaultanswer={'256','[-0.1 0.4]','gtec_import'};
        
        answer=inputdlg2(prompt,name,numlines,defaultanswer);
        
        for i=1:size(prompt,2)
            [val status] = str2num(answer{1});  % Use curly bracket for subscript
            if ~status
                % Handle empty value returned for unsuccessful conversion
                % ...
                %error('Unexpected input...');
            end
            
            % get input to variables
            switch i
                case 1
                    datasetname = answer{i};
                case 2
                    epochlimits = str2num(answer{i});
                otherwise
                    fprintf('Unhandled user input..');
                    return;
            end
        end
        
        % loading class info is documented out
        % to activate class info loading remove comment characters for the
        % whole next block!!!
        load_classinfo = questdlg2('Load event information from ... ?','Event Information','Data','File','None','Data');

        if strcmpi(load_classinfo,'File')
            [filenametmp, filepath] = uigetfile('*.mat', 'Choose a g.tec class info file -- pop_read_gtec');
            drawnow; 
        elseif strcmpi(load_classinfo,'Data')
            filenametmp = -1;
        else 
            filenametmp = 0;
        end
        
        if filenametmp ~= 0
            classinfo_filename = [filepath '/' filenametmp];   
            EEG = read_gtec_hdf5events(EEG, filename, classinfo_filename);
        elseif filenametmp == -1
            EEG = read_gtec_hdf5events(EEG, filename, 0);
        end
                
        % get sampling rate from xml
        AcquisitionTaskDescription = h5read(filename,'/RawData/AcquisitionTaskDescription');
        % get xml tree
        xmltree = xmlreadstring(AcquisitionTaskDescription);
        % get root element of xml
        root = xmltree.getDocumentElement;
        % get subelements
        elements = root.getChildNodes;
        
        % get sampling rate from xml elements
        node = elements.getFirstChild;
        while ~isempty(node)
            if strcmpi(node.getNodeName, 'SamplingFrequency')
                break;
            else
                node = node.getNextSibling;
            end
        end
        EEG.srate = str2double(node.getTextContent);

        % get number of channels
        node = elements.getFirstChild;
        while ~isempty(node)
            if strcmpi(node.getNodeName, 'NumberOfAcquiredChannels')
                break;
            else
                node = node.getNextSibling;
            end
        end
        EEG.nbchan = str2double(node.getTextContent);
        
        % get channel names
        channelNames = cell(1,EEG.nbchan);
        channelNr = 1;
        node = elements.getFirstChild;
        node = node.getNextSibling.getNextSibling.getNextSibling;
        node = node.getChildNodes;
        node = node.getFirstChild;
        while ~isempty(node)
            if strcmpi(node.getNodeName, 'ChannelProperties')
                innernode = node.getChildNodes;
                innernode = innernode.getFirstChild;
                innernode = innernode.getNextSibling;
                %innernode = innernode.getChildNodes;
                %innernode = innernode.getFirstChild;
                %innernode = innernode.getNextSibling;
                while ~isempty(innernode)
                    if strcmpi(innernode.getNodeName, 'ChannelName')
                        break;
                    else
                        innernode = innernode.getNextSibling;
                    end
                end
                channelNames(channelNr)=innernode.getTextContent;
                channelNr = channelNr + 1;
                %node = innernode.getParentNode;
                node = node.getNextSibling;
            else
                node = node.getNextSibling;
            end
        end
         
        % load channel name to channel structure
        for i=1:size(channelNames,2)
            channelLabel = channelNames{i};
            if strcmpi(channelLabel,'')
                channelLabel = num2str(i);
            end
            EEG.chanlocs(i).labels = channelLabel;
        end
        
        % import data from measurement file
        EEG.data = h5read(filename,'/RawData/Samples');
        % set dataset name
        EEG.setname = datasetname;
        
        % create epochs
        %EEG_targets = pop_epoch(EEG,num2cell(2),epochlimits,'new_name',[datasetname ' ' 'targetepochs']);
        %EEG_nontargets = pop_epoch(EEG,num2cell(1),epochlimits,'new_name',[datasetname ' ' 'nontargetepochs']);
        
        % save epochs in new datasets
        %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', datasetname);
        %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG_targets, 1, 'setname', [datasetname ' ' 'targetepochs']);
        %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG_nontargets, 1, 'setname', [datasetname ' ' 'nontargetepochs']);
        
end
a=0; 
return;
