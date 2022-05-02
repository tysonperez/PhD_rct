function varargout = mark_ICs_gui(varargin)
% MARK_ICS_GUI MATLAB code for mark_ICs_gui.fig
%      MARK_ICS_GUI, by itself, creates a new MARK_ICS_GUI or raises the existing
%      singleton*.
%
%      H = MARK_ICS_GUI returns the handle to a new MARK_ICS_GUI or the handle to
%      the existing singleton*.
%
%      MARK_ICS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARK_ICS_GUI.M with the given input arguments.
%
%      MARK_ICS_GUI('Property','Value',...) creates a new MARK_ICS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mark_ICs_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mark_ICs_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mark_ICs_gui

% Last Modified by GUIDE v2.5 13-Sep-2018 16:41:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mark_ICs_gui_OpeningFcn, ...
    'gui_OutputFcn',  @mark_ICs_gui_OutputFcn, ...
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


% --- Executes just before mark_ICs_gui is made visible.
function mark_ICs_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mark_ICs_gui (see VARARGIN)

% Choose default command line output for mark_ICs_gui
handles.output = hObject;

% set (handles.figure_mark_ICs_gui, 'toolbar', 'figure');

set(handles.button_next,'String','<html>Next &#8594</html>');
set(handles.button_previous,'String','<html>&#8592 Previous</html>');

%% input
% handles.EEG = varargin{1, :};
handles.ALLEEG = evalin ('base', 'ALLEEG');
handles.EEG = evalin ('base', 'EEG');
handles.CURRENTSET = evalin ('base', 'CURRENTSET');

handles.vertical_lines = []; % [20, 30];      % dotted vertical lines at ms latency from trigger. empty for none

handles.pvaf_threshold = 98;    % percentage to use for threshold to detect number of ICs which describes/covers at least this amount of variance in the data

%% default
handles.IC_number = 1;

handles.temporal_plot.length = 10;      % number of epochs for epoched data, or number of seconds for continous data (same protocol as EEGLAB eegplot)
handles.temporal_plot.start = 1;
handles.temporal_plot.YLim = [];
handles.temporal_plot.YLim_resolution = 1.25;    % resolution when axis is expanded/contracted. note expansion/contration in total will be 2 times of this.

handles.erp_plot.YLim = [];
handles.erp_plot.YLim_resolution = [];        % resolution when axis is expanded/contracted. note expansion/contration in total will be 2 times of this. this is dynamic, changed for every component and is set by max of abs of yLim divided by 10.

if (isfield(handles.EEG.reject.SASICA, 'icarejADJUST'))      % was ADJUST performed
    handles.ADJUST_result = handles.EEG.reject.SASICA.icarejADJUST;
else
    handles.button_adjust_unveil.Visible = 'off';
end

if (isfield(handles.EEG.reject.SASICA, 'icarejFASTER'))      % was FASTER performed
    handles.FASTER_result = handles.EEG.reject.SASICA.icarejFASTER;
else
    handles.button_faster_unveil.Visible = 'off';
end

handles.SASICA_result = handles.EEG.reject.SASICA.icarejautocorr | ...
    handles.EEG.reject.SASICA.icarejfocalcomp | ...
    handles.EEG.reject.SASICA.icarejtrialfoc | ...
    handles.EEG.reject.SASICA.icarejresvar;

if (isfield(handles.EEG.etc, 'ica_marked_component_type'))      %already marked data
    handles.component_type = handles.EEG.etc.ica_marked_component_type;
else
    for i = 1:length(handles.EEG.icachansind)
        handles.component_type(i).brain = false;
        handles.component_type(i).eye = false;
        handles.component_type(i).muscle = false;
        handles.component_type(i).heart = false;
        handles.component_type(i).line = false;
        handles.component_type(i).channel = false;
        handles.component_type(i).other = false;
        handles.component_type(i).not_sure = false;     % not_sure = false, which means i am sure.
    end
end

%MRI file
load (['dipfit_MNI', filesep, 'avg152t1.mat']);
handles.MRI = mri;

%% GUI initialization
% single dipole
handles.axes_single_dipole.sagittal = handles.axes_single_dipole_sagittal;
handles.axes_single_dipole.transverse = handles.axes_single_dipole_transverse;
handles.axes_single_dipole.coronal = handles.axes_single_dipole_coronal;

% bilateral dipole
handles.axes_bilateral_dipole.sagittal = handles.axes_bilateral_dipole_sagittal;
handles.axes_bilateral_dipole.transverse = handles.axes_bilateral_dipole_transverse;
handles.axes_bilateral_dipole.coronal = handles.axes_bilateral_dipole_coronal;

% fitTwo dipole
handles.axes_fitTwo_dipole.sagittal = handles.axes_fitTwo_dipole_sagittal;
handles.axes_fitTwo_dipole.transverse = handles.axes_fitTwo_dipole_transverse;
handles.axes_fitTwo_dipole.coronal = handles.axes_fitTwo_dipole_coronal;

% GUI subject/filter info to verify that correct file loaded
temp_filepath = strsplit(handles.EEG.filepath, ' - ');
temp_filepath = strsplit(temp_filepath{end}, filesep);
handles.filter_type = temp_filepath{1};    %FIR or IIR
handles.text_dataset_name.String = [handles.EEG.setname, ' (', handles.filter_type, ')'];

handles.editbox_temporal_plot_length.String = num2str (handles.temporal_plot.start);

handles.button_next.Enable = 'on';
handles.button_previous.Enable = 'off';

handles.show_trigger_markers = true;

handles = update_gui (handles, handles.IC_number);



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mark_ICs_gui wait for user response (see UIRESUME)
% uiwait(handles.figure_mark_ICs_gui);


% --- Outputs from this function are returned to the command line.
function varargout = mark_ICs_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_next.
function button_next_Callback(hObject, eventdata, handles)
% hObject    handle to button_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (~strcmp(handles.button_next.String, 'Save (overwrites)'))
    % at least one checkbox is marked/selected
    if (handles.checkbox_brain.Value | handles.checkbox_eye.Value | ...
            handles.checkbox_muscle.Value | handles.checkbox_heart.Value | ...
            handles.checkbox_line_noise.Value | handles.checkbox_channel_noise.Value | ...
            handles.checkbox_other.Value | handles.checkbox_not_sure.Value)
        
        if (handles.IC_number <= length(handles.EEG.icachansind))
            
            handles.IC_number = handles.IC_number + 1;
            
            handles = update_gui (handles, handles.IC_number);
            
            if (strcmp(handles.button_previous.Enable, 'off'))
                handles.button_previous.Enable = 'on';
            end
        end
    else
        errordlg({'Mark at least one checkbox'}, ...
            'Error', 'modal');
    end
    
    if (handles.IC_number == length(handles.EEG.icachansind))
        handles.button_next.String = 'Save (overwrites)';
    else
        handles.button_next.String = '<html>Next &#8594</html>';
    end
else    % save data
    button = questdlg('Save (Overwrite the file)?', 'Verify', ...
        'Yes', 'No', 'No');
    
    if (strcmp (button, 'Yes'))
        
        component_type = handles.component_type;
        
        for i = 1:length(handles.EEG.icachansind)
            is_IC_marked = struct2array(component_type(i));
            
            is_IC_marked = sum(is_IC_marked);
            
            if (is_IC_marked == 0)
                errordlg({['IC not marked: ', num2str(i)]}, ...
                    'Error', 'modal');
                
                return;
            end
        end
        
        handles.EEG.etc.ica_marked_component_type = handles.component_type;
        reject_components = false (1, length(handles.EEG.icachansind));
        
        for i = 1:length(handles.EEG.icachansind)   % reject all components which are not marked brain
            if (~handles.component_type(i).brain)
                reject_components(i) = true;
            end
        end
        handles.EEG.reject.gcompreject = reject_components;
        
        % save EEG dataset. NOTE: it overwrites the dataset.
        handles.EEG = pop_saveset(handles.EEG, 'filepath', handles.EEG.filepath, 'filename', handles.EEG.filename);
        [handles.ALLEEG, handles.EEG] = eeg_store(handles.ALLEEG, handles.EEG, handles.CURRENTSET);
        
        assignin ('base', 'ALLEEG', handles.ALLEEG);
        assignin ('base', 'EEG', handles.EEG);
        
        eeglab redraw;
        
        disp ('====================');
        disp ([handles.EEG.setname, ' : Saved']);
        disp ('====================');
        
        % plot uncleaned and cleaned data
        EEG2 = pop_subcomp(handles.EEG, find(handles.EEG.reject.gcompreject == 1), 0);
        eegplot(handles.EEG.data, 'data2', EEG2.data, 'tag', 'EEG_PLOT', 'winlength', 10, ...
            'eloc_file', handles.EEG.chanlocs, 'plottitle', [handles.EEG.filename, ' - ', 'black=uncleaned; red=cleaned']);
    end
end

guidata(hObject, handles);


% --- Executes on button press in button_previous.
function button_previous_Callback(hObject, eventdata, handles)
% hObject    handle to button_previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.IC_number > 1)
    
    handles.IC_number = handles.IC_number - 1;
    
    handles = update_gui (handles, handles.IC_number);
    
    if (strcmp(handles.button_next.String, 'Save (overwrites)'))
        handles.button_next.String = '<html>Next &#8594</html>';
    end
end

if (handles.IC_number == 1)
    handles.button_previous.Enable = 'off';
end

guidata(hObject, handles);


% --- Executes on button press in checkbox_brain.
function checkbox_brain_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_brain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_brain

if (get(hObject,'Value') == 1)
    handles.component_type(handles.IC_number).brain = true;
else
    handles.component_type(handles.IC_number).brain = false;
end

guidata(hObject, handles);


% --- Executes on button press in checkbox_eye.
function checkbox_eye_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_eye (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_eye

if (get(hObject,'Value') == 1)
    handles.component_type(handles.IC_number).eye = true;
else
    handles.component_type(handles.IC_number).eye = false;
end

guidata(hObject, handles);


% --- Executes on button press in checkbox_muscle.
function checkbox_muscle_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_muscle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_muscle

if (get(hObject,'Value') == 1)
    handles.component_type(handles.IC_number).muscle = true;
else
    handles.component_type(handles.IC_number).muscle = false;
end

guidata(hObject, handles);


% --- Executes on button press in checkbox_heart.
function checkbox_heart_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_heart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_heart

if (get(hObject,'Value') == 1)
    handles.component_type(handles.IC_number).heart = true;
else
    handles.component_type(handles.IC_number).heart = false;
end

guidata(hObject, handles);


% --- Executes on button press in checkbox_line_noise.
function checkbox_line_noise_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_line_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_line_noise

if (get(hObject,'Value') == 1)
    handles.component_type(handles.IC_number).line = true;
else
    handles.component_type(handles.IC_number).line = false;
end

guidata(hObject, handles);


% --- Executes on button press in checkbox_channel_noise.
function checkbox_channel_noise_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_channel_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_channel_noise

if (get(hObject,'Value') == 1)
    handles.component_type(handles.IC_number).channel = true;
else
    handles.component_type(handles.IC_number).channel = false;
end

guidata(hObject, handles);


% --- Executes on button press in checkbox_other.
function checkbox_other_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_other (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_other

if (get(hObject,'Value') == 1)
    handles.component_type(handles.IC_number).other = true;
else
    handles.component_type(handles.IC_number).other = false;
end

guidata(hObject, handles);


% --- Executes on button press in checkbox_not_sure.
function checkbox_not_sure_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_not_sure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_not_sure

if (get(hObject,'Value') == 1)
    handles.component_type(handles.IC_number).not_sure = true;
else
    handles.component_type(handles.IC_number).not_sure = false;
end

guidata(hObject, handles);


% --- Executes on button press in button_adjust_unveil.
function button_adjust_unveil_Callback(hObject, eventdata, handles)
% hObject    handle to button_adjust_unveil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.ADJUST_result(handles.IC_number))
    handles.button_adjust_unveil.String = 'Reject';
    handles.button_adjust_unveil.BackgroundColor = [228,26,28]/255;
else
    handles.button_adjust_unveil.String = 'Accept';
    handles.button_adjust_unveil.BackgroundColor = [77,175,74]/255;
end

handles.button_adjust_unveil.Enable = 'inactive';

guidata(hObject, handles);


% --- Executes on button press in button_faster_unveil.
function button_faster_unveil_Callback(hObject, eventdata, handles)
% hObject    handle to button_faster_unveil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.FASTER_result(handles.IC_number))
    handles.button_faster_unveil.String = 'Reject';
    handles.button_faster_unveil.BackgroundColor = [228,26,28]/255;
else
    handles.button_faster_unveil.String = 'Accept';
    handles.button_faster_unveil.BackgroundColor = [77,175,74]/255;
end

handles.button_faster_unveil.Enable = 'inactive';

guidata(hObject, handles);


% --- Executes on button press in button_sasica_unveil.
function button_sasica_unveil_Callback(hObject, eventdata, handles)
% hObject    handle to button_sasica_unveil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.SASICA_result(handles.IC_number))
    handles.button_sasica_unveil.String = 'Reject';
    handles.button_sasica_unveil.BackgroundColor = [228,26,28]/255;
else
    handles.button_sasica_unveil.String = 'Accept';
    handles.button_sasica_unveil.BackgroundColor = [77,175,74]/255;
end

handles.button_sasica_unveil.Enable = 'inactive';

guidata(hObject, handles);


% --- Executes on button press in button_next_single.
function button_next_single_Callback(hObject, eventdata, handles)
% hObject    handle to button_next_single (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

epoch_number = 1;
direction = 'inc';

handles = scroll_time_series (handles, epoch_number, direction);

guidata(hObject, handles);


% --- Executes on button press in button_next_multiple.
function button_next_multiple_Callback(hObject, eventdata, handles)
% hObject    handle to button_next_multiple (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

epoch_number = handles.temporal_plot.length;
direction = 'inc';

handles = scroll_time_series (handles, epoch_number, direction);

guidata(hObject, handles);


% --- Executes on button press in button_next_multiple_more.
function button_next_multiple_more_Callback(hObject, eventdata, handles)
% hObject    handle to button_next_multiple_more (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

epoch_number = handles.temporal_plot.length * 10;
direction = 'inc';

handles = scroll_time_series (handles, epoch_number, direction);

guidata(hObject, handles);


% --- Executes on button press in button_previous_single.
function button_previous_single_Callback(hObject, eventdata, handles)
% hObject    handle to button_previous_single (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

epoch_number = 1;
direction = 'dec';

handles = scroll_time_series (handles, epoch_number, direction);

guidata(hObject, handles);


% --- Executes on button press in button_previous_multiple.
function button_previous_multiple_Callback(hObject, eventdata, handles)
% hObject    handle to button_previous_multiple (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

epoch_number = handles.temporal_plot.length;
direction = 'dec';

handles = scroll_time_series (handles, epoch_number, direction);

guidata(hObject, handles);


% --- Executes on button press in button_previous_multiple_more.
function button_previous_multiple_more_Callback(hObject, eventdata, handles)
% hObject    handle to button_previous_multiple_more (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

epoch_number = handles.temporal_plot.length * 10;
direction = 'dec';

handles = scroll_time_series (handles, epoch_number, direction);

guidata(hObject, handles);


function editbox_temporal_plot_length_Callback(hObject, eventdata, handles)
% hObject    handle to editbox_temporal_plot_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editbox_temporal_plot_length as text
%        str2double(get(hObject,'String')) returns contents of editbox_temporal_plot_length as a double

epoch_number = str2double(get(hObject,'String'));
direction = [];

handles = scroll_time_series (handles, epoch_number, direction);

guidata(hObject, handles);


% --- Executes on button press in button_expand_amplitude.
function button_expand_amplitude_Callback(hObject, eventdata, handles)
% hObject    handle to button_expand_amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.temporal_plot.YLim = handles.temporal_plot.YLim + handles.temporal_plot.YLim_resolution;

handles = update_yAxis_time_series (handles);

guidata(hObject, handles);


% --- Executes on button press in button_compress_amplitude.
function button_compress_amplitude_Callback(hObject, eventdata, handles)
% hObject    handle to button_compress_amplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.temporal_plot.YLim = handles.temporal_plot.YLim - handles.temporal_plot.YLim_resolution;

if (handles.temporal_plot.YLim < handles.temporal_plot.YLim_resolution)
    handles.temporal_plot.YLim = handles.temporal_plot.YLim_resolution;
end

handles = update_yAxis_time_series (handles);

guidata(hObject, handles);



function editbox_temporal_plot_ylim_Callback(hObject, eventdata, handles)
% hObject    handle to editbox_temporal_plot_ylim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editbox_temporal_plot_ylim as text
%        str2double(get(hObject,'String')) returns contents of editbox_temporal_plot_ylim as a double

handles.temporal_plot.YLim = str2double(get(hObject,'String'));

if (handles.temporal_plot.YLim < handles.temporal_plot.YLim_resolution)
    handles.temporal_plot.YLim = handles.temporal_plot.YLim_resolution;
end

handles = update_yAxis_time_series (handles);

guidata(hObject, handles);


% --- Executes on button press in button_expand_amplitude_erp.
function button_expand_amplitude_erp_Callback(hObject, eventdata, handles)
% hObject    handle to button_expand_amplitude_erp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.erp_plot.YLim = handles.erp_plot.YLim + handles.erp_plot.YLim_resolution;

handles = update_yAxis_erp (handles);

guidata(hObject, handles);


% --- Executes on button press in button_compress_amplitude_erp.
function button_compress_amplitude_erp_Callback(hObject, eventdata, handles)
% hObject    handle to button_compress_amplitude_erp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.erp_plot.YLim = handles.erp_plot.YLim - handles.erp_plot.YLim_resolution;

if (handles.erp_plot.YLim < handles.erp_plot.YLim_resolution)
    handles.erp_plot.YLim = handles.erp_plot.YLim_resolution;
end

handles = update_yAxis_erp (handles);

guidata(hObject, handles);


% --- Executes on button press in button_trigger_markers.
function button_trigger_markers_Callback(hObject, eventdata, handles)
% hObject    handle to button_trigger_markers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_trigger_markers

button_state = get(hObject,'Value');
if (button_state == get(hObject,'Max'))
    handles.show_trigger_markers = false;
elseif (button_state == get(hObject,'Min'))
    handles.show_trigger_markers = true;
end

handles = plot_time_series (handles, handles.IC_number);

guidata(hObject, handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text_IC_number.
function text_IC_number_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text_IC_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

defaultans = {num2str(handles.IC_number)};

IC_number = inputdlg('Enter IC Number', '', 1, defaultans);

if (~isempty(IC_number))
    IC_number = str2double(cell2mat (IC_number));
    if ((IC_number > 0) & (IC_number <= length(handles.EEG.icachansind)) & (IC_number ~= handles.IC_number))
        handles.IC_number = IC_number;
        
        handles = update_gui (handles, handles.IC_number);
    elseif (IC_number ~= handles.IC_number)
        errordlg({'Invalid input'}, ...
            'Error', 'modal');
    end
end

guidata(hObject, handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text_IC_variance.
function text_IC_variance_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text_IC_variance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pvaf = handles.EEG.etc.ica_pvaf;

sum_ = 0;
for i = 1:length(pvaf)
   sum_ = sum_ + pvaf(i);
   
   if (sum_ >= handles.pvaf_threshold)
       break
   end
end

str_ = ['\fontsize{15}Number of ICs required for pvaf ', num2str(handles.pvaf_threshold), '% are ', num2str(i)];

CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'nonmodal';

msgbox(str_, 'ICs for pvaf threshold', 'none', CreateStruct);

guidata(hObject, handles);



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text_spectrum.
function text_spectrum_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text_spectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure;
plot (handles.EEG.etc.ica_psd_freq, handles.EEG.etc.ica_psd_spectra(handles.IC_number, :), 'LineWidth', 2);
grid on;
grid minor;
title ('Activity Power Spectrum', 'FontSize', 14);
xlabel ('Frequency (Hz)', 'FontSize', 14)
ylabel ('Log Power Spectral Density 10*log_{10}(\muV^{2}/Hz)', 'FontSize', 14);
xlim ([handles.EEG.etc.ica_psd_freq(1), handles.EEG.etc.ica_psd_freq(end)]);

guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axes_erpimage_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_erpimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nbins = handles.EEG.EVENTLIST.nbin;
bin_codes = [handles.EEG.event.bini];
bdf = {handles.EEG.EVENTLIST.bdf.description};
total_bins = [handles.EEG.EVENTLIST.eventinfo.code];

str_ = '\fontsize{15}';

for nb = 1:nbins
    block_start = find(bin_codes == nb, 1);
    block_end = find(bin_codes == nb, 1, 'last');
    block_max = find(total_bins == nb, 1, 'last') - find(total_bins == nb, 1) + 1;
    percent_accepted = 100 * (block_end - block_start + 1)/block_max;
    
    str_ = [str_, bdf{nb}, ': trials = ', num2str(block_start), ' to ', num2str(block_end), ...
        ' (', num2str(percent_accepted, 3), '%)', char(10)];
end

CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'nonmodal';

msgbox(str_, 'Block Details', 'help', CreateStruct);

guidata(hObject, handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over checkbox_not_sure.
function checkbox_not_sure_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_not_sure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

answer = questdlg('Mark the remaining ICs as?', ...
    'Mark remaining ICs', ...
    'Brain', 'Others', 'Cancel', 'Cancel');

switch answer
    case 'Brain'
        for IC_number = handles.IC_number+1:1:length(handles.EEG.icachansind)
            handles.component_type(IC_number).brain = true;
            handles.component_type(IC_number).other = false;
        end
    case 'Others'
        for IC_number = handles.IC_number+1:1:length(handles.EEG.icachansind)
            handles.component_type(IC_number).brain = false;
            handles.component_type(IC_number).other = true;
        end
    case 'Cancel'
        guidata(hObject, handles);
        
        return;
end

handles.IC_number = length(handles.EEG.icachansind);

update_gui (handles, handles.IC_number);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = update_gui (handles, IC_number)

% profile on;

% w = waitbar (0, '', 'Name', 'Processing');
% movegui (w, 'northwest');
% waitbar_steps = 1;
% waitbar_step = 0;
% waitbar_resolution = 1/8;       %1 divided by the number of times waitbar is to be updated in one loop (inner)

%% reset UI components
% waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, 'reset');
handles = reset_UI (handles, IC_number);


%% disable all clickable UI components
% waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, 'UI off');
handles = toggle_UI_components (handles, false);


%% component topoplot and component info
% waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, 'IC topoplot');
handles = plot_component (handles, IC_number);


%% erpimage
% waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, 'ERP image');
handles = plot_erpimage (handles, IC_number);


%% dipplots
% waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, 'dipplots');
handles = plot_dipplots (handles, IC_number);


%% show dipplot using EEGLAB dipplot function
% show_EEGLAB_dipplot (handles, IC_number);


%% spectrum
% the first dataset specturm was calculated at 4Hz frequency resolution
% the later ones are calculated at 1Hz. The xaxis needs to be fixed
% accordingly
% waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, 'Spectrum');
handles = plot_spectrum (handles, IC_number);


%% time series
% waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, 'Time Series');
handles = plot_time_series (handles, IC_number);


%% enable all clickable UI components
% waitbar_step = update_waitbar (w, waitbar_resolution, waitbar_step, waitbar_steps, 'UI on');
handles = toggle_UI_components (handles, true);


%%
% close (w);

% profile viewer;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = reset_UI (handles, IC_number)

%% temporal plot
handles.temporal_plot.start = 1;
handles.temporal_plot.YLim = [];

%% erp plot
handles.erp_plot.YLim = [];

%% ADJUST Unveil button
handles.button_adjust_unveil.String = 'ADJUST Unveil';
handles.button_adjust_unveil.BackgroundColor = [0.94, 0.94, 0.94];
handles.button_adjust_unveil.Enable = 'on';

%% FASTER Unveil button
handles.button_faster_unveil.String = 'FASTER Unveil';
handles.button_faster_unveil.BackgroundColor = [0.94, 0.94, 0.94];
handles.button_faster_unveil.Enable = 'on';

%% SASICA Unveil button
handles.button_sasica_unveil.String = 'SASICA Unveil';
handles.button_sasica_unveil.BackgroundColor = [0.94, 0.94, 0.94];
handles.button_sasica_unveil.Enable = 'on';

%% checkboxes
if (handles.component_type(IC_number).brain)
    handles.checkbox_brain.Value = 1;
else
    handles.checkbox_brain.Value = 0;
end

if (handles.component_type(IC_number).eye)
    handles.checkbox_eye.Value = 1;
else
    handles.checkbox_eye.Value = 0;
end

if (handles.component_type(IC_number).muscle)
    handles.checkbox_muscle.Value = 1;
else
    handles.checkbox_muscle.Value = 0;
end

if (handles.component_type(IC_number).heart)
    handles.checkbox_heart.Value = 1;
else
    handles.checkbox_heart.Value = 0;
end

if (handles.component_type(IC_number).line)
    handles.checkbox_line_noise.Value = 1;
else
    handles.checkbox_line_noise.Value = 0;
end

if (handles.component_type(IC_number).channel)
    handles.checkbox_channel_noise.Value = 1;
else
    handles.checkbox_channel_noise.Value = 0;
end

if (handles.component_type(IC_number).other)
    handles.checkbox_other.Value = 1;
else
    handles.checkbox_other.Value = 0;
end

if (handles.component_type(IC_number).not_sure)     % not_sure = false, which means i am sure.
    handles.checkbox_not_sure.Value = 1;
else
    handles.checkbox_not_sure.Value = 0;
end

%% Next/Save button check
if (IC_number == length(handles.EEG.icachansind))
    handles.button_next.String = 'Save (overwrites)';
else
    handles.button_next.String = '<html>Next &#8594</html>';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = toggle_UI_components (handles, is_enable)

if (is_enable)
    enable = 'on';
else
    enable = 'off';
end

set (findall(handles.figure_mark_ICs_gui, '-property', 'style'), 'enable', enable);

if (handles.IC_number == 1)
    handles.button_previous.Enable = 'off';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = show_EEGLAB_dipplot (handles, comp)

dipplot(handles.EEG.etc.dipfit.single_dipole.model(comp), ...
    'mri', ['dipfit_MNI', filesep, 'avg152t1.mat'], 'normlen','on', ...
    'coordformat', 'spherical', 'axistight', 'off', 'gui', 'off', 'summary', 'on');

dipplot(handles.EEG.etc.dipfit.bilateral_dipole.model(comp), ...
    'mri', ['dipfit_MNI', filesep, 'avg152t1.mat'], 'normlen','on', ...
    'coordformat', 'spherical', 'axistight', 'off', 'gui', 'off', 'summary', 'on');

dipplot(handles.EEG.etc.dipfit.two_dipole.model(comp), ...
    'mri', ['dipfit_MNI', filesep, 'avg152t1.mat'], 'normlen','on', ...
    'coordformat', 'spherical', 'axistight', 'off', 'gui', 'off', 'summary', 'on');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = plot_component (handles, IC_number)

% axes(handles.axes_topoplot);
set(handles.figure_mark_ICs_gui, 'CurrentAxes', handles.axes_topoplot)

topoplot(handles.EEG.icawinv(:,IC_number), handles.EEG.chanlocs, 'chaninfo', handles.EEG.chaninfo);

handles.text_IC_number.String = ['IC: ', num2str(IC_number), ' of ', num2str(length(handles.EEG.icachansind))];

handles.text_IC_variance.String = ['pvaf: ', num2str(handles.EEG.etc.ica_pvaf(IC_number), 3), '%'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = plot_erpimage (handles, IC_number)

vertical_lines = handles.vertical_lines;      % dotted vertical lines at ms latency from trigger

data = squeeze(handles.EEG.icaact(IC_number, :, :));

%% erp image
cla (handles.axes_erpimage);
delete(findobj('Type', 'ColorBar'))

if (IC_number == 1)
    delete(findobj('Type', 'Legend'))
end

mindat = min(min(data));
maxdat = max(max(data));
maxdat = max(abs([mindat maxdat])); % make symmetrical about 0
mindat = -maxdat;

imagesc ('Parent', handles.axes_erpimage, 'XData', handles.EEG.times, ...
    'YData', [1:1:handles.EEG.trials], 'CData', data');
handles.axes_erpimage.XLim = [handles.EEG.times(1) handles.EEG.times(end)];
handles.axes_erpimage.YLim = [1, handles.EEG.trials];
handles.axes_erpimage.CLim = [mindat, maxdat];
handles.axes_erpimage.XTick = [];
handles.axes_erpimage.LineWidth = 1.5;
handles.axes_erpimage.YAxis.LineWidth = 1.5;
handles.axes_erpimage.YAxis.TickDirection = 'out';
handles.axes_erpimage.YLabel.String = 'Trials';

% colorbar
colormap jet;
handles.axes_erp_cbar.Visible = 'off';
erp_cbar = colorbar (handles.axes_erp_cbar, 'Position', handles.axes_erp_cbar.Position);
handles.axes_erp_cbar.CLim = [mindat, maxdat];
erp_cbar.Ticks = [ceil([mindat, mindat/2]*100)/100, 0, floor([maxdat/2, maxdat]*100)/100];

% vertical lines
hold (handles.axes_erpimage, 'on');
plot (handles.axes_erpimage, [0, 0], handles.axes_erpimage.YLim, 'k', 'LineWidth', 2);  % vertical line at 0ms

if (~isempty (vertical_lines))
    for vl = 1:length(vertical_lines)
        plot (handles.axes_erpimage, [vertical_lines(vl), vertical_lines(vl)], handles.axes_erpimage.YLim, 'k--', 'LineWidth', 1.5);
    end
end
hold (handles.axes_erpimage, 'off');


% horizontal lines
nbins = handles.EEG.EVENTLIST.nbin;
bin_codes = [handles.EEG.event.bini];

hold (handles.axes_erpimage, 'on');
for hl = 1:nbins
    ydata = find(bin_codes == hl, 1);
    plot (handles.axes_erpimage, handles.axes_erpimage.XLim, [ydata, ydata], 'k-.', 'LineWidth', 0.5);
end
plot (handles.axes_erpimage, handles.axes_erpimage.XLim, [handles.axes_erpimage.YLim(2), handles.axes_erpimage.YLim(2)], 'k-.', 'LineWidth', 0.5);
hold (handles.axes_erpimage, 'off');



%% two erps
cla (handles.axes_erp_two);
handles.axes_erp_two.YLimMode = 'auto';

erp1 = data(:, 1:round(handles.EEG.trials/2));
erp2 = data(:, round(handles.EEG.trials/2)+1:end);

erp1 = mean (erp1, 2);
erp2 = mean (erp2, 2);

hold (handles.axes_erp_two, 'on');
p(1) = plot (handles.axes_erp_two, handles.EEG.times, erp1, 'Color', [228,26,28]/255, 'LineWidth', 2);
p(2) = plot (handles.axes_erp_two, handles.EEG.times, erp2, 'Color', [152,78,163]/255, 'LineWidth', 2);

yLim = max(abs(handles.axes_erp_two.YLim));
handles.axes_erp_two.YLim = [-yLim, yLim];

plot (handles.axes_erp_two, handles.axes_erp_two.XLim, [0, 0], 'k', 'LineWidth', 2);

plot (handles.axes_erp_two, [0, 0], handles.axes_erp_two.YLim, 'k', 'LineWidth', 2);
for vl = 1:length(vertical_lines)
    plot (handles.axes_erp_two, [vertical_lines(vl), vertical_lines(vl)], handles.axes_erp_two.YLim, 'k--', 'LineWidth', 1.5);
end
hold (handles.axes_erp_two, 'off');

handles.axes_erp_two.XLim = [handles.EEG.times(1) handles.EEG.times(end)];
handles.axes_erp_two.XTick = [];
handles.axes_erp_two.YLabel.String = '\muV';


%% one erp
cla (handles.axes_erp_one);
handles.axes_erp_one.YLimMode = 'auto';

erp1 = mean (data, 2);

hold (handles.axes_erp_one, 'on');
p(3) = plot (handles.axes_erp_one, handles.EEG.times, erp1, 'Color', [0, 153, 26]/255, 'LineWidth', 2);

% yLim = max(abs(handles.axes_erp_one.YLim));
handles.axes_erp_one.YLim = [-yLim, yLim];

plot (handles.axes_erp_one, handles.axes_erp_one.XLim, [0, 0], 'k', 'LineWidth', 2);

plot (handles.axes_erp_one, [0, 0], handles.axes_erp_one.YLim, 'k', 'LineWidth', 2);
for vl = 1:length(vertical_lines)
    plot (handles.axes_erp_one, [vertical_lines(vl), vertical_lines(vl)], handles.axes_erp_one.YLim, 'k--', 'LineWidth', 1.5);
end
hold (handles.axes_erp_one, 'off');

handles.axes_erp_one.XLim = [handles.EEG.times(1) handles.EEG.times(end)];
handles.axes_erp_one.YLabel.String = '\muV';
handles.axes_erp_one.XLabel.String = 'Time (ms)';

%% y-axis
% update y-axis
if (isempty (handles.erp_plot.YLim))
    handles.erp_plot.YLim = max (abs(handles.axes_erp_two.YLim));
    handles.erp_plot.YLim_resolution = max (abs(handles.axes_erp_two.YLim))/10;
end
handles = update_yAxis_erp (handles);

%%legend
if (IC_number == 1)
    legend (p, {'Lower 50%', 'Upper 50%', 'All Trials'}, 'FontSize', 10, 'Location', 'northwest');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = plot_dipplots (handles, IC_number)

% reset axes for dipplots before using them again
field_planes = {'sagittal', 'transverse', 'coronal'};
for idx = 1:length(field_planes)
    s = getfield(handles.axes_single_dipole, field_planes{idx});
    b = getfield(handles.axes_bilateral_dipole, field_planes{idx});
    f = getfield(handles.axes_fitTwo_dipole, field_planes{idx});
    
    cla(s, 'reset');
    cla(b, 'reset');
    cla(f, 'reset');
end

% single dipole
ica_dipplot(handles.EEG.etc.dipfit.single_dipole.model(IC_number), handles.axes_single_dipole, ...
    'mri', handles.MRI, 'normlen','on', ...
    'coordformat', 'spherical', 'axistight', 'off', 'gui', 'off', 'summary', 'on', 'verbose', 'off');

handles.text_single_dipole_RV.String = ['RV: ', num2str(100*handles.EEG.etc.dipfit.single_dipole.model(IC_number).rv, 3), '%'];


% bilateral dipole
ica_dipplot(handles.EEG.etc.dipfit.bilateral_dipole.model(IC_number), handles.axes_bilateral_dipole, ...
    'mri', handles.MRI, 'normlen','on', ...
    'coordformat', 'spherical', 'axistight', 'off', 'gui', 'off', 'summary', 'on', 'verbose', 'off');

handles.text_bilateral_dipole_RV.String = ['RV: ', num2str(100*handles.EEG.etc.dipfit.bilateral_dipole.model(IC_number).rv, 3), '%'];


% fitTwo dipole
if (handles.EEG.etc.dipfit.fitTwo_dipole_success && isfield(handles.EEG.etc.dipfit.fitTwo_dipole.model, 'select'))
    if (isempty (handles.EEG.etc.dipfit.fitTwo_dipole.model(IC_number).select))
        handles.text_fitTwo_dipole.Visible = 'off';
        handles.panel_fitTwo_dipole.Visible = 'off';
        handles.text_fitTwo_dipole_RV.Visible = 'off';
    else
        handles.text_fitTwo_dipole.Visible = 'on';
        handles.panel_fitTwo_dipole.Visible = 'on';
        handles.text_fitTwo_dipole_RV.Visible = 'on';
        
        ica_dipplot(handles.EEG.etc.dipfit.fitTwo_dipole.model(IC_number), handles.axes_fitTwo_dipole, ...
            'mri', handles.MRI, 'normlen','on', ...
            'coordformat', 'spherical', 'axistight', 'off', 'gui', 'off', 'summary', 'on', 'verbose', 'off');
        
        handles.text_fitTwo_dipole_RV.String = ['RV: ', num2str(100*handles.EEG.etc.dipfit.fitTwo_dipole.model(IC_number).rv, 3), '%'];
    end
else
    handles.text_fitTwo_dipole.Visible = 'off';
    handles.panel_fitTwo_dipole.Visible = 'off';
    handles.text_fitTwo_dipole_RV.Visible = 'off';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = plot_spectrum (handles, IC_number)

% the first dataset, healthy FIR, spectrum was calculated at 4Hz frequency
% resolution. The later ones are to be done on 1Hz. The xaxis needs to be
% fixed accordingly.
% for 4Hz resolution: 0 to 40 Hz (ica_psd_freq array, the indices for 0 to
% 40 Hz = 1:11)
% for 1Hz resolution: 0 to 40 Hz (ica_psd_freq array, the indices for 0 to
% 40 Hz = 1:41)
plot (handles.axes_spectrum1, handles.EEG.etc.ica_psd_freq(1:41),handles.EEG.etc.ica_psd_spectra(IC_number, 1:41), 'r', 'LineWidth', 2)
handles.axes_spectrum1.XTick = [0:10:40];
grid (handles.axes_spectrum1, 'on')
set (handles.axes_spectrum1, 'FontSize', 14);
handles.axes_spectrum1.YLabel.String = '10*log_{10}(\muV^{2}/Hz)';

% for 4Hz resolution: 0 to 80 Hz (ica_psd_freq array, the indices for 0 to
% 80 Hz = 1:21)
% for 1Hz resolution: 0 to 80 Hz (ica_psd_freq array, the indices for 0 to
% 80 Hz = 1:81)
plot (handles.axes_spectrum2, handles.EEG.etc.ica_psd_freq(1:81), handles.EEG.etc.ica_psd_spectra(IC_number, 1:81), 'r', 'LineWidth', 2)
handles.axes_spectrum2.XTick = [0:10:80];
grid (handles.axes_spectrum2, 'on')
set (handles.axes_spectrum2, 'FontSize', 14);
handles.axes_spectrum2.XLabel.String = 'Frequency (Hz)';
handles.axes_spectrum2.YLabel.String = 'Log Power Spectral Density_.';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = plot_time_series (handles, IC_number)

cla (handles.axes_temporal);
handles.axes_temporal.YLimMode = 'auto';

hold (handles.axes_temporal, 'on');

if (~isempty(handles.EEG.epoch))    % data is epoched
    
    trials_to_plot = handles.temporal_plot.start:1:handles.temporal_plot.start+handles.temporal_plot.length-1;
    
    plt = squeeze (handles.EEG.icaact (IC_number, :, trials_to_plot));
    plt = reshape (plt, [1, size(plt,1)*size(plt,2)]);      % make continous signal by concatenating epochs
    
    plot (handles.axes_temporal, plt);
    
    handles.axes_temporal.XLim = [1, length(plt)];
    handles.axes_temporal.XTick = 1:handles.EEG.pnts/2:length(plt);
    handles.axes_temporal.XTick = handles.axes_temporal.XTick(2:2:end);
    handles.axes_temporal.XTickLabel = trials_to_plot;
    
    trigger_markers = find(handles.EEG.times == 0);     % trigger is at 0ms
    trigger_markers = [trigger_markers:handles.EEG.pnts:length(plt)];
    epoch_markers = trigger_markers - trigger_markers(1) + 1;       % epoch start/end markers / this separates epochs
    
    % epoch start/end marker
    plot (handles.axes_temporal, [epoch_markers; epoch_markers], ...
        [handles.axes_temporal.YLim(1); handles.axes_temporal.YLim(2)], '--', 'Color', [99, 99, 99]/255);
    
    % trigger marker
    if (handles.show_trigger_markers)
        plot (handles.axes_temporal, [trigger_markers; trigger_markers], ...
            [handles.axes_temporal.YLim(1); handles.axes_temporal.YLim(2)], 'r');
    end
    
    % update y-axis
    if (isempty (handles.temporal_plot.YLim))
        handles.temporal_plot.YLim = max (abs(handles.axes_temporal.YLim));
    end
    handles = update_yAxis_time_series (handles);
    
    handles.temporal_plot.start = trials_to_plot(1);
    
    handles.editbox_temporal_plot_length.String = num2str (handles.temporal_plot.start);
    
else    % data is continuous
    
end

hold (handles.axes_temporal, 'off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = scroll_time_series (handles, epoch_number, direction)

if (isempty(direction))     % make the left most epoch number on plot equal to the editbox input number, (if it falls in ploting range)
    epoch_number = epoch_number;
elseif (strcmp (direction, 'inc'))
    epoch_number = handles.temporal_plot.start + epoch_number;
elseif (strcmp (direction, 'dec'))
    epoch_number = handles.temporal_plot.start - epoch_number;
end

if (epoch_number <= 0)
    epoch_number = 1;
elseif (epoch_number > (handles.EEG.trials - handles.temporal_plot.length + 1))     % can't go any furhter towards right
    epoch_number = handles.EEG.trials - handles.temporal_plot.length + 1;
end

handles.temporal_plot.start = epoch_number;

handles = plot_time_series (handles, handles.IC_number);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = update_yAxis_time_series (handles)

handles.axes_temporal.YLim = [-handles.temporal_plot.YLim, handles.temporal_plot.YLim];

vertical_lines = handles.axes_temporal.Children;
for i = 1:1:length(vertical_lines)-1
    vertical_lines(i).YData =[handles.axes_temporal.YLim(1); handles.axes_temporal.YLim(2)];
end

handles.editbox_temporal_plot_ylim.String = num2str (handles.temporal_plot.YLim, 4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = update_yAxis_erp (handles)

handles.axes_erp_two.YLim = [-handles.erp_plot.YLim, handles.erp_plot.YLim];
handles.axes_erp_one.YLim = [-handles.erp_plot.YLim, handles.erp_plot.YLim];

vertical_lines_erp_two = handles.axes_erp_two.Children;
vertical_lines_erp_one = handles.axes_erp_one.Children;
for i = 1:1:length(vertical_lines_erp_two)-3
    vertical_lines_erp_two(i).YData =[handles.axes_erp_two.YLim(1); handles.axes_erp_two.YLim(2)];
    vertical_lines_erp_one(i).YData =[handles.axes_erp_one.YLim(1); handles.axes_erp_one.YLim(2)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
