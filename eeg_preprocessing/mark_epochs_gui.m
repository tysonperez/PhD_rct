function varargout = mark_epochs_gui(varargin)
% MARK_EPOCHS_GUI MATLAB code for mark_epochs_gui.fig
%      MARK_EPOCHS_GUI, by itself, creates a new MARK_EPOCHS_GUI or raises the existing
%      singleton*.
%
%      H = MARK_EPOCHS_GUI returns the handle to a new MARK_EPOCHS_GUI or the handle to
%      the existing singleton*.
%
%      MARK_EPOCHS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARK_EPOCHS_GUI.M with the given input arguments.
%
%      MARK_EPOCHS_GUI('Property','Value',...) creates a new MARK_EPOCHS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mark_epochs_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mark_epochs_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mark_epochs_gui

% Last Modified by GUIDE v2.5 02-Jun-2018 03:46:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mark_epochs_gui_OpeningFcn, ...
    'gui_OutputFcn',  @mark_epochs_gui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mark_epochs_gui is made visible.
function mark_epochs_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mark_epochs_gui (see VARARGIN)

% Choose default command line output for mark_epochs_gui
handles.output = hObject;

%% input
% handles.EEG = varargin{1, :};
handles.ALLEEG = evalin ('base', 'ALLEEG');
handles.EEG = evalin ('base', 'EEG');
handles.CURRENTSET = evalin ('base', 'CURRENTSET');

%% bits from ERPLAB and main_preprocessing
% variables contain the bit number which is 1 of the binary.
% bits are ordered 8 to 1 (instead of 7 to 0). Bit 1 is LSB, and always
% set/high if the epoch is marked artifact. Bits 8 to 2 define the type(s)
% of artifacts. (See ERPLAB help for bit details)
bit.artifacts = 1;                  % 0000 0001
bit.max_voltage_threshold = 2;      % 0000 0010
bit.pp_threshold = 3;               % 0000 0100
bit.step_like = 4;                  % 0000 1000
bit.sample_difference = 5;          % 0001 0000
bit.flatline = 6;                   % 0010 0000
bit.others = 8;                     % 1000 0000

bit.combination = bit.artifacts;    % will contain the artifact types, using bitor to enter them and bitand to remove.

handles.bit = bit;

%% EVENTLIST arrays
handles.EVENTLIST.original = cell2mat({handles.EEG.EVENTLIST.eventinfo.flag});
handles.EVENTLIST.modified = cell2mat({handles.EEG.EVENTLIST.eventinfo.flag});


%% marker arrays
handles.markers.original = handles.EEG.reject.rejmanual;
handles.markers.modified = handles.EEG.reject.rejmanual;



%% default values

% GUI subject/filter info to verify that correct file loaded
temp_filepath = strsplit(handles.EEG.filepath, ' - ');
temp_filepath = strsplit(temp_filepath{end}, filesep);
handles.filter_type = temp_filepath{1};    %FIR or IIR
set (handles.text_dataset_name, 'String', [handles.EEG.setname, ' (', handles.filter_type, ')']);

handles.marker_type = 'modified';       % original or modified;
set (handles.radio_modified, 'Value', 1);
handles.EVENTLIST.selected = handles.EVENTLIST.modified;

handles.boolean_operator = 'OR';       % AND or OR;
set (handles.radio_OR, 'Value', 1);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mark_epochs_gui wait for user response (see UIRESUME)
% uiwait(handles.figure_mark_epochs_gui);


% --- Outputs from this function are returned to the command line.
function varargout = mark_epochs_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on button press in checkbox_maxVoltageThreshold.
function checkbox_maxVoltageThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_maxVoltageThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_maxVoltageThreshold

if (get(hObject,'Value') == get(hObject,'Max'))
    tmp = bitset (1, handles.bit.max_voltage_threshold, 1);
    handles.bit.combination = bitor (handles.bit.combination, tmp);
else
    tmp = bitset (255, handles.bit.max_voltage_threshold, 0);
    handles.bit.combination = bitand (handles.bit.combination, tmp);
end

guidata(hObject, handles);


% --- Executes on button press in checkbox_ppVoltageThreshold.
function checkbox_ppVoltageThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ppVoltageThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ppVoltageThreshold
if (get(hObject,'Value') == get(hObject,'Max'))
    tmp = bitset (1, handles.bit.pp_threshold, 1);
    handles.bit.combination = bitor (handles.bit.combination, tmp);
else
    tmp = bitset (255, handles.bit.pp_threshold, 0);
    handles.bit.combination = bitand (handles.bit.combination, tmp);
end

guidata(hObject, handles);


% --- Executes on button press in checkbox_step.
function checkbox_step_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_step
if (get(hObject,'Value') == get(hObject,'Max'))
    tmp = bitset (1, handles.bit.step_like, 1);
    handles.bit.combination = bitor (handles.bit.combination, tmp);
else
    tmp = bitset (255, handles.bit.step_like, 0);
    handles.bit.combination = bitand (handles.bit.combination, tmp);
end

guidata(hObject, handles);


% --- Executes on button press in checkbox_sampleDifference.
function checkbox_sampleDifference_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_sampleDifference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_sampleDifference
if (get(hObject,'Value') == get(hObject,'Max'))
    tmp = bitset (1, handles.bit.sample_difference, 1);
    handles.bit.combination = bitor (handles.bit.combination, tmp);
else
    tmp = bitset (255, handles.bit.sample_difference, 0);
    handles.bit.combination = bitand (handles.bit.combination, tmp);
end

guidata(hObject, handles);


% --- Executes on button press in checkbox_flatline.
function checkbox_flatline_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_flatline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_flatline
if (get(hObject,'Value') == get(hObject,'Max'))
    tmp = bitset (1, handles.bit.flatline, 1);
    handles.bit.combination = bitor (handles.bit.combination, tmp);
else
    tmp = bitset (255, handles.bit.flatline, 0);
    handles.bit.combination = bitand (handles.bit.combination, tmp);
end

guidata(hObject, handles);


% --- Executes on button press in checkbox_others.
function checkbox_others_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_others (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_others
if (get(hObject,'Value') == get(hObject,'Max'))
    tmp = bitset (1, handles.bit.others, 1);
    handles.bit.combination = bitor (handles.bit.combination, tmp);
else
    tmp = bitset (255, handles.bit.others, 0);
    handles.bit.combination = bitand (handles.bit.combination, tmp);
end

guidata(hObject, handles);


% --- Executes when selected object is changed in radiogroup_markers.
function radiogroup_markers_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in radiogroup_markers
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get (eventdata.NewValue, 'Tag')
    case 'radio_original'
        handles.marker_type = 'original';
        handles.EVENTLIST.selected = handles.EVENTLIST.original;
        
    case 'radio_modified'
        handles.marker_type = 'modified';
        handles.EVENTLIST.selected = handles.EVENTLIST.modified;
end

guidata(hObject, handles);


% --- Executes when selected object is changed in radiogroup_boolean_operator.
function radiogroup_boolean_operator_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in radiogroup_boolean_operator
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get (eventdata.NewValue, 'Tag')
    case 'radio_AND'
        handles.boolean_operator = 'AND';
        
    case 'radio_OR'
        handles.boolean_operator = 'OR';
end

guidata(hObject, handles);


% --- Executes on button press in pushbutton_plot.
function pushbutton_plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_9plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.bit.combination <= 1)
    errordlg('Select at least one type of artifact', 'Error');
else
    handles = toggle_UI_components (handles, false);
    
    if (get(handles.radio_modified, 'Value') == 1)      % check here needed to get updated modified eventlist if the radiobutton is not alternated to and fro after running plot for modified.
        handles.EVENTLIST.selected = handles.EVENTLIST.modified;
    else
        handles.EVENTLIST.selected = handles.EVENTLIST.original;
    end
    
    % get the epochs marked rejected, of the types of artifacts selected
    % (ANDed or ORed) from the selected (original or modified) marker array
    idx = getArtifactIndices (handles.EVENTLIST.selected, handles.bit.combination, handles.boolean_operator);
    handles.EEG.reject.rejmanual = idx;
    
    % change channel names for easy identification of interpolated and bad
    % channels given by PREP
    chanlocs = handles.EEG.chanlocs;
    interpolatedChannels = handles.EEG.etc.noiseDetection.interpolatedChannelNumbers;
    stillNoisyChannels = handles.EEG.etc.noiseDetection.stillNoisyChannelNumbers;
    
    for i = 1:length(interpolatedChannels)
        handles.EEG.chanlocs(interpolatedChannels(i)).labels = ['i--', handles.EEG.chanlocs(interpolatedChannels(i)).labels];
    end
    
    for i = 1:length(stillNoisyChannels)
        handles.EEG.chanlocs(stillNoisyChannels(i)).labels = ['n--', handles.EEG.chanlocs(stillNoisyChannels(i)).labels];
    end
    
    % plot
    assignin ('base', 'EEG', handles.EEG);      % overwrite EEG struct in base workspace, as it will be overwritten by the eegplot. Overwriting it here will overwrite the EEG.reject.rejmanual (which is equal to idx), and only that will be processed in setArtifactIndices.
    
    pop_eegplot(handles.EEG, 1, 1, 0, ...
        0, 'tag', 'EEG_PLOT', ...
        'winlength', 10, ...
        'color', {[0,0,0]/255, [55,126,184]/255, [77,175,74]/255, [152,78,163]/255, [255,127,0]/255});  % give colors (4-8) from color brewer site?
    
    while (~isempty(findobj('tag','EEG_PLOT')))
        pause (5);      % changed from 0.5 to 2, after seeing good scrolling of comparison EEG sets after ICA cleaning
    end
    
    % always update the marker modified array, keeping the original array
    % untouched
    handles.EEG = evalin ('base', 'EEG');   %since changes are written in the base workspace by the eegplot, get the EEG struct.
    
    [handles.EVENTLIST.modified, handles.markers.modified] = ...
        setArtifactIndices (handles.EVENTLIST.modified, idx, ...
        handles.bit.combination, handles.bit.others, handles.EEG.reject);
    
    % put back original in base workspace
    handles.EEG.reject.rejmanual = handles.markers.original;     % put back the original markers
    handles.EEG.chanlocs = chanlocs;            % original channel info struct (channel labels)
    assignin ('base', 'EEG', handles.EEG);      % now EEG struct in base workspace will be the original EEG struct.
    
    handles = toggle_UI_components (handles, true);
end

guidata(hObject, handles);


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

button = questdlg('Save (Overwrite the file)?', 'Verify', ...
    'Yes', 'No', 'No');

if (strcmp (button, 'Yes'))
    
    handles = toggle_UI_components (handles, false);
    
    % commented for the reason mentioned in setArtifactIndices
    % bit_others = bitset (1, handles.bit.others, 1);
    %
    % idx_manual = getArtifactIndices (handles.EVENTLIST.modified, bit_others, 'OR');
    % idx_artifacts = getArtifactIndices (handles.EVENTLIST.modified, 255, 'OR');
    %
    % idx = idx_manual | idx_artifacts;
    %
    % handles.EEG.reject.rejmanualE = repmat (idx, size(handles.EEG.reject.rejmanualE));
    
    % get subject number
    filename = strsplit(handles.EEG.filename, '-');
    filename = filename{1};
    
    % 15-apr-2019
    % to have both the subject number and session number in the file
    % instead of just the session number (represented orignally by the
    % filename)
    temp_subject = strsplit(handles.EEG.filepath, filesep);
    temp_subject = temp_subject{end-3};
    temp_subject = strsplit(temp_subject, ' ');
    subject_number = str2double(temp_subject{2});
    
    % update EVENTLIST flag struct
    for i = 1:length(handles.EVENTLIST.modified)
        handles.EEG.EVENTLIST.eventinfo(i).flag = handles.EVENTLIST.modified(i);
        handles.EEG.epoch(i).eventflag = handles.EVENTLIST.modified(i);     % added 06-march-2018. EEG.epoch is not updated when erplab2eeglab is used to syncronize, so manually did here.
    end
    
    % transfer EVENTLIST to REJECT using erplab.
    handles.EEG = pop_syncroartifacts(handles.EEG, 'Direction', 'erplab2eeglab');
    
    % save the event list txt for modified markers
    handles.EEG = pop_exporteegeventlist(handles.EEG , ...
        'Filename', [handles.EEG.filepath, filename, '-eventList-ARTmanual.txt']);
    
    
    % save stats in mat file
    main_folder = strsplit(handles.EEG.filepath, filesep);
    main_folder = main_folder (1:end-4);
    main_folder = strjoin (main_folder, filesep);
    
    [EEG, tprej, acce, rej, histoflags] = pop_summary_AR_eeg_detection(handles.EEG, [handles.EEG.filepath, filename, '-stats-ARTmanual.txt']);   %prints table to file
    [~, ~, ~, ~, ~] = pop_summary_AR_eeg_detection(handles.EEG, '');   %prints table on command line  % added 06-march-2018 since separate_manually_marked_artifacts script is commented below as EEG.epoch sync is fixed manually above.
    
    load ([main_folder, filesep, 'data_info', filesep, 'epoch_artifacts_stats_', handles.filter_type, '.mat']);
    
    if (~exist('stats_ARTmanual', 'var'))
        stats_ARTmanual = [];
    end
    
    %     ALLEEG = handles.ALLEEG;                  % commented 06-march-2018
    %     EEG = handles.EEG;                        % commented 06-march-2018
    %     separate_manually_marked_artifacts;       % commented 06-march-2018
    number_rejected_trials = length(find(handles.EVENTLIST.modified ~= 0));
    
    stats_ARTmanual(length(stats_ARTmanual)+1).subject = subject_number;
    stats_ARTmanual(length(stats_ARTmanual)).session = str2double(filename);
    stats_ARTmanual(length(stats_ARTmanual)).accepted = length(handles.EVENTLIST.modified) - number_rejected_trials;
    stats_ARTmanual(length(stats_ARTmanual)).rejected = number_rejected_trials;  %total epochs rejected after visual inspection
    stats_ARTmanual(length(stats_ARTmanual)).percent_rejected = 100*number_rejected_trials/length(handles.EVENTLIST.modified);
    stats_ARTmanual(length(stats_ARTmanual)).F8_others = histoflags(:, 8);
    clear ALLEEG EEG;
    
    save ([main_folder, filesep, 'data_info', filesep, 'epoch_artifacts_stats_', handles.filter_type, '.mat'], 'stats_ARTmanual', '-append');
    
    % get difference of original and modified markers and put in EEG etc struct
    handles.EEG.etc.EVENTLIST.modified = handles.EVENTLIST.modified;
    handles.EEG.etc.EVENTLIST.reject_difference = xor (handles.EVENTLIST.original, handles.EVENTLIST.modified);     % difference b/w original and modified markers
    
    % save EEG dataset. NOTE: it overwrites the dataset.
    handles.EEG = pop_saveset(handles.EEG, 'savemode','resave');
    [handles.ALLEEG, handles.EEG] = eeg_store(handles.ALLEEG, handles.EEG, handles.CURRENTSET);
    
    assignin ('base', 'ALLEEG', handles.ALLEEG);
    assignin ('base', 'EEG', handles.EEG);
    
    eeglab redraw;
    
    disp ('====================');
    disp ([handles.EEG.setname, ' : Saved']);
    disp ('====================');
    
    handles = toggle_UI_components (handles, true);
end

guidata(hObject, handles);


% --- Executes when user attempts to close figure_mark_epochs_gui.
function figure_mark_epochs_gui_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure_mark_epochs_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if (strcmp(get(handles.pushbutton_save, 'enable'), 'on'))
    delete(hObject);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%
function [idx] = getArtifactIndices (eventinfo, bits_artifact, boolean_operator)
% input:    eventinfo: EEG.EVENTLIST.eventinfo
%           bits_artifact: Decimal of the binary combination
%           (e.g. 5 for 0000 0101, 7 for 0000 0111)
% output:   idx: logical index array containing 1s for the type of artifacts
%           given as input to be plotted

idx = bitand (eventinfo, bits_artifact);

if (strcmp (boolean_operator, 'AND'))
    idx = (idx == bits_artifact);
elseif (strcmp (boolean_operator, 'OR'))    % separate each type of artifact
    tmp = dec2bin(bits_artifact, 8);
    
    tmp = cellstr (tmp');
    
    tmp = str2double (tmp)';
    
    tmp = fliplr (tmp);
    tmp = find (tmp == 1);
    
    idx2 = false(size (idx));
    
    for i = 2:1:length(tmp)
        tmp(i) = flag2dec ('ArtifactFlag', [1, tmp(i)]);
        
        tmp_idx = bitand (eventinfo, tmp(i));
        tmp_idx = (tmp_idx == tmp(i));
        
        idx2 = idx2 | tmp_idx;
    end
    
    idx = idx2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% the idea was to have bit_others (bit 8 in the first version of code) set
% to 1 and bit_artifact (bit 1 in the first version of code) to 0 (as it is
% 1 if the epoch is marked artifact by ERPLAB) to identify that manually
% unmarked the artifact epoch. For marked artifact, set bit_others (bit 8)
% to 1 and bit_artifact (bit 1) to 1, to identify manual marking of
% artifact. However, the unmarking strategy didn't work, as ERPLAB sets
% the epoch as artifact if any bit is 1 when ERPLAB eventlist and EEG
% reject struct are synched in the save function. So, now just putting 0 in
% the unmarked artifact. It can be used to see the difference b/w original
% and manually modified markers when running this gui. The difference array
% is saved in the EEG etc struct as well for later need if required.

function [eventinfo, markers] = setArtifactIndices (eventinfo, idx, ...
    bits_artifact, bits_others, EEG_reject)

idx_changed = xor (idx, EEG_reject.rejmanual);  % get epochs changed i.e manually marked and unmarked

artifacts_unmarked = idx & idx_changed;  % epochs manually unmarked for removal from automated marking (i.e. marked as correct epochs)

artifacts_marked = xor (idx_changed, artifacts_unmarked);    % epochs manually marked as artifacts

bits_others_unmarked = bitset (0, bits_others, 1);
% tmp = bitxor (eventinfo(artifacts_unmarked == 1), bits_artifact);
% eventinfo (artifacts_unmarked == 1) = bitor (tmp, bits_others_unmarked);
eventinfo (artifacts_unmarked == 1) = 0;

bits_others_marked = bitor (bits_others_unmarked, 1);
tmp = bitxor (eventinfo(artifacts_marked == 1), bits_artifact);
eventinfo (artifacts_marked == 1) = bitor (tmp, bits_others_marked);

markers = EEG_reject.rejmanual;


%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = toggle_UI_components (handles, is_enable)

if (is_enable)
    enable = 'on';
else
    enable = 'off';
end

set (findall(handles.figure_mark_epochs_gui, '-property', 'style'), 'enable', enable);

handles.UI_enable = is_enable;

