% pop_loadhdf5() - Read g-tec 24bit EEG file from g-recorder
%
% Usage:
%   >> EEG = pop_loadhdf5;                                                 % an interactive uigetfile window
%   >> EEG = pop_loadhdf5( 'filename','filename.hdf5', filepath, 'c:/');   % no pop-up window
%   >> EEG = pop_loadhdf5( PropertyName, PropertyValue);
%
% Optional properties:
%   filename       - Filename of the *.hdf5 file
%   filepath       - File path of the *-hdf5 file
%   rejectchans    - Channel numbers to reject (remove) from the data
%   ref_ch         - specify the reference channel or channels. (type: string cell array, size: Nref x 1)
%                    This should be the channel that was plugged in to the reference input of the amplifier - the data will not be referenced!
%   ref_range      - Ranges of channels to assign the ref_ch.
%                    The number of ranges must corrospond to the number of reference channels. (type: doubble cell array, size: Nref x 1)
%                    Eg. pop_loadhdf5('ref_ch',{'TP9','TP7','Cz'},'ref_range',{4:7,8:11,12:43})
%
% Outputs:
%   EEG            - EEGLAB data structure
%
% Author: Simon Lind Kappel, Aarhus University, 2015
%
% Revision 1.10  2015/04/20 - Simon
%  Fix issues with uicontrol
%  Add properties to define reference channels.
%  Add uicontrol field to specify the reference channel (same reference for all channels)


function [EEG, command] = pop_loadhdf5(varargin)

command = '';
EEG = [];

%% Create input argument structure
if ~isempty(varargin)
    try
        g = struct(varargin{:});
    catch
        disp('Setevent: wrong syntax in function arguments');
        return;
    end
else
    g = struct();
end;

%% History
if ~isempty(which('vararg2str'))
    command = sprintf('[EEG, command] = pop_loadhdf5(%s);', vararg2str(varargin));
end

%% Create new EEG struct
if ~isempty(which('eeg_emptyset'))
    EEG = eeg_emptyset;
end

%% test the presence of variables
if isfield(g, 'filepath')
    EEG.filepath = g(1).filepath;
else
    EEG.filepath = pwd;
end

if ~isfield(g, 'filename')
    [EEG.filename, EEG.filepath, FilterIndex] = uigetfile({'*.hdf5','*.hdf5 (gRecorder)'},'Select a recording from the g-tec amplifier');
    if FilterIndex < 1, return, end
else
    EEG.filename = g(1).filename;
end

if (~strcmp(EEG.filename(end-4:end),'.hdf5'))
    EEG.filename = strcat(EEG.filename,'.hdf5');
end

Filename = fullfile(EEG.filepath,EEG.filename);

%% Read the hdf5 file and extract the EEG data and other options

%Read info about all attributes in the hdf5 file.
fileinfo = hdf5info(Filename);

%Read Raw EEG data and scale from µV to V
[EEG.data, ~] = hdf5read(Filename,'/RawData/Samples');
EEG.data = single(EEG.data);

%Get the XML data in the AcquisitionTaskDescription attribute
[AcquisitionTaskDescription, ~] = hdf5read(Filename,'/RawData/AcquisitionTaskDescription');

%Extract samplerate
FsTemp = ReadXmlAttribute(AcquisitionTaskDescription.Data, 'SampleRate');
EEG.srate = str2double(FsTemp{1});

%Extract Channel number
ChNrTemp = ReadXmlAttribute(AcquisitionTaskDescription.Data, 'PhysicalChannelNumber');
%EEG.gtec.chnumber = str2double(ChNrTemp);

%Extract Channel number
DeviceNumberTemp = ReadXmlAttribute(AcquisitionTaskDescription.Data, 'DeviceNumber');
%EEG.gtec.devnumber = str2double(DeviceNumberTemp);

%Get names of the EEG channels
ChannelNames = ReadXmlAttribute(AcquisitionTaskDescription.Data, 'ChannelName');

%Find the array index of the '/AsynchronData' item.
AsynchronDataIdx = 0;
for n=1:length(fileinfo.GroupHierarchy.Groups)
    if strcmp(fileinfo.GroupHierarchy.Groups(n).Name,'/AsynchronData')
        AsynchronDataIdx = n;
    end
end

%% Extract the trigger data from the hdf5 file and create the event struct for EEGlab

%The trigger data is stored in '/AsynchronData' [ Groups(1) ]. Check if any data is
%stores in this directory.
if (AsynchronDataIdx > 0)
    if (length(fileinfo.GroupHierarchy.Groups(AsynchronDataIdx).Datasets) > 1)
        
        %Get the XML data in the '/AsynchronData/AsynchronSignalTypes' attribute
        [TriggerChannelInfo, ~] = hdf5read(Filename,'/AsynchronData/AsynchronSignalTypes');
        TriggerName = ReadXmlAttribute(TriggerChannelInfo.Data, 'Name');
        
        %TrigChannelNum = ReadXmlAttribute(TriggerChannelInfo.Data, 'ChannelNumber');
        TrigID = str2double(ReadXmlAttribute(TriggerChannelInfo.Data, 'ID'));
        
        [TrigTime, ~] = hdf5read(Filename,'/AsynchronData/Time');
        [TrigTypeID, ~] = hdf5read(Filename,'/AsynchronData/TypeID');
        
        if length(TrigTime) ~= length(TrigTypeID)
            warning(sprintf('Trigtype is missing for the last %i triggers',length(TrigTime)-length(TrigTypeID)));
            TrigTime = TrigTime(1:length(TrigTypeID));
        end
        
        EEG.event = struct('type',{},'position',[],'latency',[],'urevent',[]); %Columns: 1=type, 2=position, 3=latency, 4=urevent
        [~,TrigSort] = sort(TrigTime);
        
        TrigTime = double(TrigTime);
        
        for n=1:length(TrigTime);
            EEG.event(n).type = TriggerName{TrigID == TrigTypeID(TrigSort(n))};
            EEG.event(n).position = 1;
            EEG.event(n).latency = TrigTime(TrigSort(n));
            EEG.event(n).urevent = n;
        end
    end
end

%% handle input arguments
rejectchanss = [];
ref_ch = [];
ref_range = [];
tmpfields = [];

if ~isempty(g)
    tmpfields = fieldnames(g);
end

for curfield = tmpfields'
    switch lower(curfield{1})
        case {'filename' } % do nothing now
        case {'filepath' } % do nothing now
        case {'rejectchans'}
            rejectchanss = getfield(g, {1}, curfield{1});
        case {'ref_ch'}
            ref_ch = {g.(curfield{1})};
        case {'ref_range'}
            ref_range = {g.(curfield{1})};
        otherwise, error(['pop_editset() error: unrecognized field ''' curfield{1} '''']);
    end;
end;

if ~any(ismember(lower(tmpfields),'rejectchans')) && ~any(ismember(lower(tmpfields),'ref_ch'))
    if isempty(which('inputgui')) == false
        res = inputgui( 'geometry', [2 2], ...
            'geomvert', [1 1], 'uilist', { ...
            { 'style', 'text', 'string', [ 'Reject channels:' ] }, ...
            { 'style', 'edit', 'string', num2str(rejectchanss)},...
            { 'style', 'text', 'string', [ 'Reference channel:' ] }, ...
            { 'style', 'edit', 'string', num2str(ref_ch)}});
        if ~isempty(res{1})
            try
                rejectchanss = eval(res{1});
            catch
                warning('The reject channel string was not formatted correctly');
            end
        end
        if ~isempty(res{2})
            ref_ch = res{2};
            if ~iscell(ref_ch)
                ref_ch = {ref_ch};
            end
        end
    else
        warning('The ''inputgui()'' EEGlab function was not located. The GUI to specify rejection channels and reference could not be displayed.')
    end
end

EEG.chanlocs = struct('labels', cellstr(ChannelNames));
EEG.ref = 'unknown';
if ~isempty(ref_ch)
    for n = 1:length(EEG.chanlocs)
        if isempty(ref_range) && length(ref_ch) == 1
            if iscell(ref_ch)
                EEG.chanlocs(n).ref = ref_ch{1};
                EEG.ref = ref_ch{1};
            elseif ischar(ref_ch)
                EEG.chanlocs(n).ref = ref_ch;
            end
        elseif length(ref_range) == length(ref_ch) && iscell(ref_range) && iscell(ref_ch)
            for refidx = 1:length(ref_range)
                if (ismember(n,ref_range{refidx}) == true)
                    EEG.chanlocs(n).ref = ref_ch{refidx};
                end
            end
        else
            warning('The ref_range or ref_ch was not formatted correctly. They must be cell arrays.')
        end
    end
end

if ~isempty(rejectchanss) && isempty(find(isnan(rejectchanss) == 1,1))
    range_chan = 1:size(EEG.data,1);
    range_chan(rejectchanss) = [];
    EEG.data = EEG.data(range_chan,:);
    EEG.chanlocs = EEG.chanlocs(range_chan);
    %EEG.gtec.chnumber = EEG.gtec.chnumber(range_chan);
    %EEG.gtec.devnumber = EEG.gtec.devnumber(range_chan);
end

%% Asign field in the EEG struct
EEG.nbchan          = size(EEG.data,1);
EEG.pnts            = size(EEG.data,2);
EEG.trials          = 1;
EEG.xmin            = 0;
EEG.xmax            = EEG.xmin + (EEG.pnts-1)/EEG.srate;
EEG.setname 		= EEG.filename(1:end-5);

if ~isempty(which('eeg_checkset'))
    EEG = eeg_checkset(EEG);
end

%% Extract data from Xml string
function [AttributeValue] = ReadXmlAttribute(DataString, AttributeName)
AttributeLocated = 0;
AttributeValueIndex = 1;
temp = '';
for i=length(AttributeName)+3:length(DataString)-length(AttributeName)-3
    
    % Start of string has been located
    if strcmp(DataString(i-length(AttributeName)-2:i-1), ['<' AttributeName '>'])
        AttributeLocated = 1;
        temp = '';
    end
    
    % End of string has been located
    if strcmp(DataString(i:i+length(AttributeName)+2), ['</' AttributeName '>'])
        AttributeLocated = 0;
        AttributeValue{AttributeValueIndex} = temp;
        AttributeValueIndex = AttributeValueIndex + 1;
    end
    
    if (AttributeLocated)
        temp(AttributeLocated) = DataString(i);
        AttributeLocated = AttributeLocated +1;
    end
end

if (~exist('AttributeValue','var'))
    AttributeValue{1} = '';
end
