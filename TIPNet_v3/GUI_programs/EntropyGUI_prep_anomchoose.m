function varargout = EntropyGUI_prep_anomchoose(varargin)
% ENTROPYGUI_PREP_ANOMCHOOSE MATLAB code for EntropyGUI_prep_anomchoose.fig
%      ENTROPYGUI_PREP_ANOMCHOOSE, by itself, creates a new ENTROPYGUI_PREP_ANOMCHOOSE or raises the existing
%      singleton*.
%
%      H = ENTROPYGUI_PREP_ANOMCHOOSE returns the handle to a new ENTROPYGUI_PREP_ANOMCHOOSE or the handle to
%      the existing singleton*.
%
%      ENTROPYGUI_PREP_ANOMCHOOSE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENTROPYGUI_PREP_ANOMCHOOSE.M with the given input arguments.
%
%      ENTROPYGUI_PREP_ANOMCHOOSE('Property','Value',...) creates a new ENTROPYGUI_PREP_ANOMCHOOSE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EntropyGUI_prep_anomchoose_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EntropyGUI_prep_anomchoose_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EntropyGUI_prep_anomchoose

% Last Modified by GUIDE v2.5 20-Jan-2016 10:04:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EntropyGUI_prep_anomchoose_OpeningFcn, ...
                   'gui_OutputFcn',  @EntropyGUI_prep_anomchoose_OutputFcn, ...
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

% --- Executes just before EntropyGUI_prep_anomchoose is made visible.
function EntropyGUI_prep_anomchoose_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EntropyGUI_prep_anomchoose (see VARARGIN)

% Choose default command line output for EntropyGUI_prep_anomchoose
handles.output = hObject;

load('GUI_programs/Temps/selection.mat');
Selected = Selected-1;
Selected

load('GUI_programs/Temps/projectname.mat')
load(fullpathname) %project file: contains mi structure of model information, data_orig, data

set(handles.enter_anom,'String',num2str(mi.DataPrep.anom.anom));
set(handles.enter_dt,'String',num2str(mi.DataPrep.anom.dt));
if strcmp(mi.DataPrep.anom.dt_unit,'mins')==1
    set(handles.button_dt,'selectedobject',handles.choose_dtmins);
elseif strcmp(mi.DataPrep.anom.dt_unit,'hours')==1
    set(handles.button_dt,'selectedobject',handles.choose_dthours);
elseif strcmp(mi.DataPrep.anom.dt_unit,'days')==1
    set(handles.button_dt,'selectedobject',handles.choose_dtdays);
elseif strcmp(mi.DataPrep.anom.dt_unit,'months')==1
    set(handles.button_dt,'selectedobject',handles.choose_dtmonths);
end
 
if strcmp(mi.DataPrep.anom.anom_unit,'days')==1
set(handles.button_anom,'selectedobject',handles.choose_anomdays);
elseif strcmp(mi.DataPrep.anom.anom_unit,'years')==1
set(handles.button_anom,'selectedobject',handles.choose_anomyears);   
end

% %update dataset according to normalization and outlier choices
% load('GUI_programs/Temps/selection.mat');
% 
% load('GUI_programs/Temps/projectname.mat')
% load(fullpathname) %project file: contains mi structure of model information, data_orig, data

size(mi.DataPrep.dp)

data = mi.DataPrep.dp{1,Selected}.data_orig;
axes(handles.plot_orig)
plot(data)

%take anomaly according to presets
mi.DataPrep.dp{Selected}.data_process1 = AnomalyFunction(mi.DataPrep.dp{Selected}.data_orig,...
    mi.DataPrep.anom.dt, mi.DataPrep.anom.dt_unit,...
    mi.DataPrep.anom.anom,mi.DataPrep.anom.anom_unit);

mi.DataPrep.dp{Selected}.data_process2 = mi.DataPrep.dp{Selected}.data_process1;

axes(handles.plot_anom)
plot(mi.DataPrep.dp{Selected}.data_process2)

handles.mi = mi;

% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = EntropyGUI_prep_anomchoose_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function enter_anom_Callback(hObject, eventdata, handles)

mi = handles.mi;
anom = str2double(get(hObject,'String'));

mi.DataPrep.anom.anom = anom; %update structure

handles.mi = mi;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function enter_anom_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function enter_dt_Callback(hObject, eventdata, handles)

mi = handles.mi;
dt = str2double(get(hObject,'String'));
mi.DataPrep.anom.dt = dt; %update structure

handles.mi = mi;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function enter_dt_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in button_dt.
function button_dt_SelectionChangeFcn(hObject, eventdata, handles)

mi = handles.mi;
Choice=get(eventdata.NewValue,'Tag'); % Get Tag of selected object

switch Choice % Get Tag of selected object.
    case 'choose_dtmins'
        dt_units = 'mins';
    case 'choose_dthours'
        dt_units = 'hrs';
    case 'choose_dtdays'
        dt_units = 'days';
    case 'choose_dtmonths'
        dt_units = 'months';
end

mi.DataPrep.anom.dt_unit = dt_units; %update structure
handles.mi = mi;

guidata(hObject,handles)


% --- Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)

%reset processed data to original data
load('GUI_programs/Temps/selection.mat');
Selected = Selected -1;
Selected

load('GUI_programs/Temps/projectname.mat')

load(fullpathname)

mi.DataPrep.dp{Selected}.data_process1 = mi.DataPrep.dp{Selected}.data_orig;
mi.DataPrep.dp{Selected}.data_process2 = mi.DataPrep.dp{Selected}.data_orig;
mi.DataPrep.alldata_process1(:,Selected) = mi.DataPrep.alldata_orig(:,Selected);
mi.DataPrep.alldata_process2(:,Selected) = mi.DataPrep.alldata_orig(:,Selected);

save(fullpathname,'mi','-append')

activefig = findobj(0,'Tag','EntropyGUI_prep');
panel = findobj(0,'Tag','select_filtering');
choice =findobj(0,'Tag','button_none');
set(panel,'selectedobject',choice);
set(activefig,'Visible','on');


close EntropyGUI_prep_anomchoose

% --- Executes on button press in button_done.
function button_done_Callback(hObject, eventdata, handles)

%update dataset according to normalization and outlier choices
load('GUI_programs/Temps/selection.mat');
Selected = Selected -1;
Selected
load('GUI_programs/Temps/projectname.mat')

mi = handles.mi;


newdata = AnomalyFunction(mi.DataPrep.dp{Selected}.data_orig,...
    mi.DataPrep.anom.dt, mi.DataPrep.anom.dt_unit,...
    mi.DataPrep.anom.anom,mi.DataPrep.anom.anom_unit);

mi.DataPrep.dp{Selected}.data_process2 = newdata;
mi.DataPrep.dp{Selected}.data_process1 = newdata;

save(fullpathname,'mi','-append')

close EntropyGUI_prep_anomchoose
activefig = findobj(0,'Tag','EntropyGUI_prep');
set(activefig,'Visible','On');


% --- Executes when selected object is changed in button_anom.
function button_anom_SelectionChangeFcn(hObject, eventdata, handles)

Choice=get(eventdata.NewValue,'Tag'); % Get Tag of selected object

mi = handles.mi;

switch Choice % Get Tag of selected object.
    case 'choose_anomdays'
        anom_units = 'days';
    case 'choose_anomyears'
        anom_units = 'years';
end

load('GUI_programs/Temps/selection.mat');

load('GUI_programs/Temps/projectname.mat')
load(fullpathname) %project file: contains mi structure of model information, data_orig, data

mi.DataPrep.anom.anom_unit = anom_units; %update structure

handles.mi = mi;
guidata(hObject,handles)


% --- Executes on button press in compute_anom.
function compute_anom_Callback(hObject, eventdata, handles)

mi = handles.mi;

load('GUI_programs/Temps/selection.mat');
Selected = Selected -1

mi.DataPrep.dp{Selected}.data_process1 = AnomalyFunction(mi.DataPrep.dp{Selected}.data_orig,...
    mi.DataPrep.anom.dt, mi.DataPrep.anom.dt_unit,mi.DataPrep.anom.anom,mi.DataPrep.anom.anom_unit);

axes(handles.plot_anom)
plot(mi.DataPrep.dp{Selected}.data_process1)

handles.mi = mi;
guidata(hObject,handles)
