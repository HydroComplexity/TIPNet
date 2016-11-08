function varargout = EntropyGUI_mainwindow(varargin)
% ENTROPYGUI_MAINWINDOW MATLAB code for EntropyGUI_mainwindow.fig
% Last Modified by GUIDE v2.5 08-Nov-2016 14:06:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EntropyGUI_mainwindow_OpeningFcn, ...
                   'gui_OutputFcn',  @EntropyGUI_mainwindow_OutputFcn, ...
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

% --- Executes just before EntropyGUI_mainwindow is made visible.
function EntropyGUI_mainwindow_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EntropyGUI_mainwindow (see VARARGIN)

clc
setappdata(hObject, 'IgnoreCloseAll', 1)
close all

addpath(genpath('Functions'));
addpath(genpath('GUI_programs'));

if exist('GUI_programs/Temps/projectname.mat','file')>0
delete('GUI_programs/Temps/projectname.mat')
end

if exist('GUI_programs/Temps/gendata.mat','file')>0
delete('GUI_programs/Temps/gendata.mat')
end

if exist('GUI_programs/Temps/dataset.mat','file')>0
delete('GUI_programs/Temps/dataset.mat')
end

%disable options until file loaded
set(handles.DataPreprocess,'Enable','Off')
set(handles.SetParams,'Enable','Off')
set(handles.set_pdfoptions,'Enable','Off')
set(handles.RunEntropyCode,'Enable','Off')
set(handles.PlotResults,'Enable','Off')

% Choose default command line output for EntropyGUI_mainwindow
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EntropyGUI_mainwindow wait for user response (see UIRESUME)
% uiwait(handles.EntropyGUI_main);

% --- Outputs from this function are returned to the command line.
function varargout = EntropyGUI_mainwindow_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on button press in LoadProjectFile.
function LoadProjectFile_Callback(hObject, eventdata, handles)

[filename1, filepath]=uigetfile({'*.mat'},...
  'Select Project File');

if filename1 ~= 0
d = load([filepath filename1]);
end

if isfield(d,'data')==1 && isfield(d,'varnames')==1 
    if isfield(d,'mi')==1
    namestring = filename1;
    fullpathname = [filepath,namestring];
    save('GUI_programs/Temps/projectname.mat','namestring','fullpathname')
    msgbox('project file loaded');
    
    %make figure handles visible for processing and network comps
    set(handles.DataPreprocess,'Enable','On')
    set(handles.SetParams,'Enable','On')
    set(handles.set_pdfoptions,'Enable','On')
    set(handles.RunEntropyCode,'Enable','On')
    
    if isfield(d,'entropy')==1
    set(handles.PlotResults,'Enable','On')
    else
    set(handles.PlotResults,'Enable','Off')      
    end
    
    else
        msgbox('not a project file, load as new dataset','error','error')
    end
         
else
    msgbox('invalid project file','error','error')
end

guidata(hObject,handles)


% --- Executes on button press in PlotResults.
function PlotResults_Callback(hObject, eventdata, handles)

load('GUI_programs/Temps/projectname.mat')
load(fullpathname) %project file: contains mi structure of model information, data_orig, data


if exist('AllStats')==0
    h=msgbox('need to compute links first','error','error');
else
EntropyGUI_plot;
wind = findobj(0,'Tag','EntropyGUI_main');
set(wind,'Visible','Off')
end


% --- Executes on button press in RunEntropyCode.
function RunEntropyCode_Callback(hObject, eventdata, handles)

load('GUI_programs/Temps/projectname.mat')
load(fullpathname)

mi.Range = [min(mi.DataPrep.alldata_process2); max(mi.DataPrep.alldata_process2)];
save(fullpathname,'mi','-append');
%Run network code
    %addpath('D:\Users\Allison\Entropy\Functions');
    
    if exist('entropy','var')
    clear entropy
    end
    
    %save(fullpathname,'entropy','-append')
    if mi.parallel_opt == 0 || mi.nsegs < 2 %no parallel
        
        for i=1:mi.nsegs
            netdata{i}.data_network = netdata{i}.data_process2;
            dat = netdata{i}.data_network;
            entropy{i}=EntropyFun_18GUI(mi,dat,i);
            save(fullpathname,'entropy','-append');     
        end
        
    elseif mi.parallel_opt == 1
    
        parfor_progress(mi.nsegs);
%         for i =1:mi.nsegs
%             entropy{i}=1;
%         end
        parfor i=1:mi.nsegs
            dat = netdata{i}.data_process2;
            entropy{i}=EntropyFun_18GUI(mi,dat,i);
            parfor_progress;
        end
        
    parfor_progress(0);
    save(fullpathname,'entropy','-append');
    end
    
    AllStats = TotalStatsFunctionGUI(netdata,mi,entropy);
    save(fullpathname,'AllStats','-append')

    h=msgbox(['TIPNet Computations Completed, Results saved as ', fullpathname]);
      
set(handles.PlotResults,'Enable','On')
    
% --- Executes on button press in SetParams.
function SetParams_Callback(hObject, eventdata, handles)

EntropyGUI_params;
wind = findobj(0,'Tag','EntropyGUI_main');
set(wind,'Visible','Off')


% --- Executes on button press in LoadNewData.
function LoadNewData_Callback(hObject, eventdata, handles)

EntropyGUI_loadnewdata
wind = findobj(0,'Tag','EntropyGUI_main');
set(wind,'Visible','Off')
waitfor(wind,'Visible')

%enable or disable buttons depending on if data was loaded
if exist('GUI_programs/Temps/dataset.mat','file')>0 &&...
        exist('GUI_programs/Temps/projectname.mat','file')>0
set(handles.DataPreprocess,'Enable','On')
set(handles.SetParams,'Enable','On')
set(handles.set_pdfoptions,'Enable','On')
set(handles.RunEntropyCode,'Enable','On')
else
set(handles.DataPreprocess,'Enable','Off')
set(handles.SetParams,'Enable','Off')
set(handles.set_pdfoptions,'Enable','Off')
set(handles.RunEntropyCode,'Enable','Off')  
end

% --- Executes on button press in GenerateData.
function GenerateData_Callback(hObject, eventdata, handles)
% open window EntropyGUI_GenerateData to create test dataset
EntropyGUI_gendata
wind = findobj(0,'Tag','EntropyGUI_main');
set(wind,'Visible','Off')

% --- Executes on button press in DataPreprocess.
function DataPreprocess_Callback(hObject, eventdata, handles)
% open window for filtering and segmenting datasets
EntropyGUI_preprocess

wind = findobj(0,'Tag','EntropyGUI_main');
set(wind,'Visible','Off')


% --- Executes on button press in set_pdfoptions.
function set_pdfoptions_Callback(hObject, eventdata, handles)
% hObject    handle to set_pdfoptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
EntropyGUI_pdfoptions

wind = findobj(0,'Tag','EntropyGUI_main');
set(wind,'Visible','Off')


% --- Executes on button press in button_exit.
function button_exit_Callback(hObject, eventdata, handles)
close(gcf)
