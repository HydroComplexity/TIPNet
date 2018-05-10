function varargout = EntropyGUI_loadnewdata(varargin)
% ENTROPYGUI_LOADNEWDATA MATLAB code for EntropyGUI_loadnewdata.fig
%      ENTROPYGUI_LOADNEWDATA, by itself, creates a new ENTROPYGUI_LOADNEWDATA or raises the existing
%      singleton*.
%
%      H = ENTROPYGUI_LOADNEWDATA returns the handle to a new ENTROPYGUI_LOADNEWDATA or the handle to
%      the existing singleton*.
%
%      ENTROPYGUI_LOADNEWDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENTROPYGUI_LOADNEWDATA.M with the given input arguments.
%
%      ENTROPYGUI_LOADNEWDATA('Property','Value',...) creates a new ENTROPYGUI_LOADNEWDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EntropyGUI_loadnewdata_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EntropyGUI_loadnewdata_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EntropyGUI_loadnewdata

% Last Modified by GUIDE v2.5 21-Jan-2016 11:10:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EntropyGUI_loadnewdata_OpeningFcn, ...
                   'gui_OutputFcn',  @EntropyGUI_loadnewdata_OutputFcn, ...
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


% --- Executes just before EntropyGUI_loadnewdata is made visible.
function EntropyGUI_loadnewdata_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

cla(handles.plot_newdata,'reset')

set(hObject,'toolbar','figure')
set(handles.button_done,'Visible','Off')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EntropyGUI_loadnewdata wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = EntropyGUI_loadnewdata_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


% --- Executes on button press in load_newmatfile.
function load_newmatfile_Callback(hObject, eventdata, handles)
%Load in a matfile, look for field 'data' as minimum requirement
%also may have fields 'timestep' and 'varnames'

clear filename1
clear data
clear d

[filename1 pathname]=uigetfile({'*.mat'},...
  'Select Data File');

if filename1 ~= 0
    d = load([pathname filename1]);


if isfield(d,'data')==1             %data is a field
    data = d.data;
    
    if isfield(d,'varnames')==1     %varnames are included
        varnames = d.varnames;
        
        if isfield(d,'timestep')==1 %timestep is included
            timestep = d.timestep;
        else
            timestep = 1:size(data,1);
        end
    else                            %varnames not included
        
        varnames = num2cell(1:size(data,2));
        
        if isfield(d,'timestep')==1 %timestep is included
            timestep = d.timestep;
        else
            timestep = 1:size(data,1);
        end
             
    end
    
nvars = size(data,2);
nSteps = size(data,1);

cla(handles.plot_newdata,'reset')
axes(handles.plot_newdata)
hold on
if nSteps < 10000
plot(timestep,data)
xlabel('time steps')
else 
plot(timestep(1:10000),data(1:10000,:))
xlabel(sprintf('time steps (first 10,000 shown of %d',nSteps))
end
legend(varnames,'Location','EastOutside')
title('Data Loaded')
xlabel('time steps')

save('GUI_programs/Temps/dataset.mat','varnames','data','nvars','nSteps','timestep')
set(handles.button_done,'Visible','On')

else
msgbox('invalid input file, must have "data" variable','error','error');
end

else %no data entered
    msgbox('no data entered','warning','warn');
end

guidata(hObject,handles)

% --- Executes on button press in load_excelfile.
function load_excelfile_Callback(hObject, eventdata, handles)


[filename1,filepath1]=uigetfile({'*.x*',},...
  'Select Excel Data File');

disp(filepath1)

addpath(filepath1)

if filename1 ~= 0
[data,varnames,raw] = xlsread(filename1);

timestep = data(:,1);
data(:,1)=[];

nvars = size(data,2);
nSteps = size(data,1);
vnames = [];

for i=1:nvars
   vnames{i}=varnames{1,i+1}; 
end

varnames=vnames;

cla(handles.plot_newdata,'reset')
axes(handles.plot_newdata)
hold on
plot(timestep,data)
legend(varnames,'Location','EastOutside')
title('Data Loaded')
xlabel('time')

save('GUI_programs/Temps/dataset.mat','varnames','data','nvars','nSteps','timestep')

set(handles.button_done,'Visible','On')

else %no data entered
    msgbox('no data entered','warning','warn');       
end

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function plot_newdata_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in button_done.
function button_done_Callback(hObject, eventdata, handles)

%save new project file and set parameters to default values
load('GUI_programs/Temps/dataset.mat')

%set all parameters and options as defaults
data(isnan(data))=0;
data_orig = data; %original dataset (data may be later altered)

[FileName,PathName] = uiputfile('*.mat','Save Project As:');

namestring = FileName;
fullpathname = [PathName FileName];

save('GUI_programs/Temps/projectname.mat','namestring','fullpathname')

%%%%%%% set all model information (mi) values to default values %%%%%%
mi.nvars = nvars;
mi.nSteps = nSteps;
mi.segsteps = nSteps; %initialize as one segment
mi.nsegs = 1;
mi.netopt = 3; %default to run full network (1=H only, 2=I,H only)
mi.nTests = 100;
mi.N = 25;
mi.nlags = 10;
mi.lagvect = 1:mi.nlags;
mi.bin_scheme = 'global';
mi.maxS = mi.nvars*3;
mi.Range = [min(data); max(data)];
mi.method = 'fixed';
mi.parallel_opt = 0;
mi.ZeroLagOpt = 0;
mi.NoSelfOpt = 0;


mi.DataPrep.FilterOpt = ones(1,nvars+1);
mi.DataPrep.NormOpt = zeros(1,nvars+1);
mi.DataPrep.OutlierOpt = zeros(1,nvars+1);
mi.DataPrep.FilterType = ones(1,nvars+1);
mi.DataPrep.FilterLambda = round(ones(1,nvars+1).*mi.nSteps./50);
mi.DataPrep.Z_effect = zeros(1,nvars+1);

for i =1:nvars
   mi.DataPrep.dp{i}.data_orig=data(:,i);
   mi.DataPrep.dp{i}.data_process1=data(:,i);
   mi.DataPrep.dp{i}.data_process2=data(:,i);
end

mi.DataPrep.alldata_process1 = data;
mi.DataPrep.alldata_process2 = data;
mi.DataPrep.alldata_orig = data;

netdata{1}.data_orig = data;
netdata{1}.data_process1 = data;
netdata{1}.data_process2 = data;
netdata{1}.timestep = timestep;

mi.DataPrep.anom.anom=5;
mi.DataPrep.anom.dt=60;
mi.DataPrep.anom.anom_unit = 'days';
mi.DataPrep.anom.dt_unit = 'mins';




guidata(hObject,handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(fullpathname,'data','mi','data_orig','varnames','netdata','timestep');

close EntropyGUI_loadnewdata
wind = findobj(0,'Tag','EntropyGUI_main');
set(wind,'Visible','On')



% --- Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)

%clear any created datasets
if exist('GUI_programs/Temps/dataset.mat')>0
delete('GUI_programs/Temps/dataset.mat')
end

close EntropyGUI_loadnewdata

wind = findobj(0,'Tag','EntropyGUI_main');
set(wind,'Visible','On')
