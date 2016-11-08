function varargout = EntropyGUI_gendata(varargin)
% ENTROPYGUI_GENDATA MATLAB code for EntropyGUI_gendata.fig
%      ENTROPYGUI_GENDATA, by itself, creates a new ENTROPYGUI_GENDATA or raises the existing
%      singleton*.
%
%      H = ENTROPYGUI_GENDATA returns the handle to a new ENTROPYGUI_GENDATA or the handle to
%      the existing singleton*.
%
%      ENTROPYGUI_GENDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENTROPYGUI_GENDATA.M with the given input arguments.
%
%      ENTROPYGUI_GENDATA('Property','Value',...) creates a new ENTROPYGUI_GENDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EntropyGUI_gendata_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EntropyGUI_gendata_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EntropyGUI_gendata

% Last Modified by GUIDE v2.5 22-Jan-2016 10:57:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EntropyGUI_gendata_OpeningFcn, ...
                   'gui_OutputFcn',  @EntropyGUI_gendata_OutputFcn, ...
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


% --- Executes just before EntropyGUI_gendata is made visible.
function EntropyGUI_gendata_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.

% Choose default command line output for EntropyGUI_gendata
handles.output = hObject;

set(hObject,'toolbar','figure')

%initialize with 2 node case, chaotic logistic feedback, full coupling
%strength and no noise
set(handles.text_eps,'String',num2str(1));
set(handles.text_epsz,'String',num2str(0));

set(handles.gendata_type,'selectedobject',handles.generate_2chaotic_ind);

set(handles.slider_eps,'Min',0);
set(handles.slider_eps,'Max',1);
set(handles.slider_eps,'Value',1);
set(handles.slider_epsz,'Min',0);
set(handles.slider_epsz,'Max',1);
set(handles.slider_epsz,'Value',0);

set(handles.generate_inputlags,'Enable','Off')

%generate chaotic logistic data
mi.N=25;
mi.eps=1;
mi.epsz=0;
mi.nvars = 2;
mi.nSteps = 1000;
varnames ={'X1','X2'};

mi.gentype = 1; %type of generated data (1-5 choices)
mi.lagmat = [0 2; 5 0]; %feedback linked

set(handles.enter_nsteps,'String',num2str(mi.nSteps));

if exist('GUI_programs/Temps/gendata.mat')>0
delete('GUI_programs/Temps/gendata.mat')
end

%generate original dataset
data = chaotic_generate(mi.nSteps,mi.lagmat,mi.eps,mi.epsz);
data_orig=data;

save('GUI_programs/Temps/gendata.mat','mi','data','data_orig','varnames')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EntropyGUI_gendata wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EntropyGUI_gendata_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_epsz_Callback(hObject, eventdata, handles)

h = get(hObject,'Value');
set(handles.text_epsz,'String',num2str(h))

load('GUI_programs/Temps/gendata.mat')

mi.epsz=h;

save('GUI_programs/Temps/gendata.mat','mi','-append')

% --- Executes during object creation, after setting all properties.
function slider_epsz_CreateFcn(hObject, eventdata, handles)

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_eps_Callback(hObject, eventdata, handles)

h = get(hObject,'Value');
set(handles.text_eps,'String',num2str(h))

load('GUI_programs/Temps/gendata.mat')

mi.eps=h;

save('GUI_programs/Temps/gendata.mat','mi','-append')


% --- Executes during object creation, after setting all properties.
function slider_eps_CreateFcn(hObject, eventdata, handles)

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in gendata.
function gendata_Callback(hObject, eventdata, handles)

load('GUI_programs/Temps/gendata.mat')
nSteps = mi.nSteps;
opt = mi.gentype;
eps = mi.eps;
epsz = mi.epsz;
lagmat = mi.lagmat;
mi.nvars=2;
mi.netopt=3; %full network default

if opt ==1 %feedback
    lagmat = [0 2; 5 0];
    Data = chaotic_generate(nSteps,lagmat,eps,epsz);
elseif opt==2 %driven (X1-->X1,X2)
    lagmat = [1 0; 5 0];
    Data = chaotic_generate(nSteps,lagmat,eps,epsz);
elseif opt ==3 %random driver
    Data(:,1:2) = rand(nSteps,2); 
    for j=3:nSteps
        Data(j,2)=(1-eps)*4*Data(j-1,2)*(1-Data(j-1,2))+...
            eps*(1-epsz)*4*Data(j-2,1)*(1-Data(j-2,1))+...
            epsz*rand();
    end
    
elseif opt==4  %read from input lag matrix
    
    Data = chaotic_generate(nSteps,lagmat,eps,epsz);
    mi.nvars = size(lagmat,1);
end

cla(handles.plot_gendata,'reset')
axes(handles.plot_gendata)
plot(Data);

nvars = mi.nvars;
data = Data;
data_orig = Data;

for i =1:mi.nvars
    varnames{i}=['X_',num2str(i)];
end

save('GUI_programs/Temps/gendata.mat','mi','data','data_orig','varnames','-append')


% --- Executes on button press in button_save.
function button_save_Callback(hObject, eventdata, handles)

load('GUI_programs/Temps/gendata.mat')

[FileName,PathName] = uiputfile('*.mat');
namestring = FileName;
fullpathname = [PathName FileName];

save('GUI_programs/Temps/projectname.mat','namestring','fullpathname');

%%%% set all model information values as defaults %%%%%%
mi.nTests = 100;
mi.N = 25;
mi.nlags = 10;
mi.lagvect = 1:mi.nlags;
mi.bin_scheme = 'global';
mi.maxS = mi.nvars*3;
mi.Range = [min(data); max(data)];
mi.segsteps = mi.nSteps;
mi.nsegs=1;
mi.method = 'KDE';
mi.parallel_opt = 0;

mi.DataPrep.FilterOpt = ones(1,mi.nvars+1);
mi.DataPrep.NormOpt = zeros(1,mi.nvars+1);
mi.DataPrep.OutlierOpt = zeros(1,mi.nvars+1);

mi.DataPrep.dp{1}.data_orig = data;
mi.DataPrep.dp{1}.data_process1 = data;
mi.DataPrep.dp{1}.data_process2 = data;

for i =1:mi.nvars
   mi.DataPrep.dp{i+1}.data_orig=data(:,i);
   mi.DataPrep.dp{i+1}.data_process1=data(:,i);
   mi.DataPrep.dp{i+1}.data_process2=data(:,i);
end

netdata{1}.data_orig = mi.DataPrep.dp{1}.data_orig;
netdata{1}.data_process1 = mi.DataPrep.dp{1}.data_process1;
netdata{1}.data_process2 = mi.DataPrep.dp{1}.data_process2;

mi.DataPrep.anom.anom=5;
mi.DataPrep.anom.dt=60;
mi.DataPrep.anom.anom_unit = 'days';
mi.DataPrep.anom.dt_unit = 'mins';

data_orig=data;

mi.DataPrep.dp{1}.data_orig = data;
mi.DataPrep.dp{1}.data_process1 = data;
mi.DataPrep.dp{1}.data_process2 = data;

for i=1:mi.nvars
  
  %varterm: variance of non-zero values
  ndata = sum(mi.DataPrep.dp{i+1}.data_process2~=0);
  
  mi.KDEparams.h1(i)= 1.06 * ndata^(-1/5) .* var(mi.DataPrep.dp{i+1}.data_process2); 
  mi.KDEparams.h2(i)= 1.77 * ndata^(-1/6) .* var(mi.DataPrep.dp{i+1}.data_process2); 
  mi.KDEparams.h3(i) = 2.78 * ndata^(-1/7) .* var(mi.DataPrep.dp{i+1}.data_process2); 
end
mi.KDEparams.h1_seg_orig = mi.KDEparams.h1; 
mi.KDEparams.h2_seg_orig=mi.KDEparams.h2;
mi.KDEparams.h3_seg_orig=mi.KDEparams.h3;

save(fullpathname,'data','mi','data_orig','varnames','netdata');

close EntropyGUI_gendata
wind = findobj(0,'Tag','EntropyGUI_main');
set(wind,'Visible','On')

% --- Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)


close EntropyGUI_gendata
wind = findobj(0,'Tag','EntropyGUI_main');
set(wind,'Visible','On')

% --- Executes when selected object is changed in gendata_type.
function gendata_type_SelectionChangeFcn(hObject, eventdata, handles)

Choice=get(eventdata.NewValue,'Tag'); % Get Tag of selected object
load('GUI_programs/Temps/gendata.mat')

switch Choice % Get Tag of selected object.
    case 'generate_2chaotic_feedback'
        mi.gentype = 1;
    case 'generate_2chaotic_driven'
        mi.gentype = 2;
    case 'generate_2chaotic_randomX'
        mi.gentype = 3;
    case 'generate_inputlags'
        mi.gentype = 4;
        
        [filename1]=uigetfile({'*.mat'},...
            'Select Lagged Adjacency Matrix');
        
        if filename1 ~= 0
            d = load(filename1);
        else
            mi.gentype =1;
        end
        
        if isfield(d,'lagmat')==1
            mi.lagmat = d.lagmat;
        end
end

if mi.gentype==1
    set(hObject,'Tag','generate_2chaotic_feedback')
end

save('GUI_programs/Temps/gendata.mat','mi','-append')


function enter_nsteps_Callback(hObject, eventdata, handles)

nSteps = str2double(get(hObject,'String'));
load('GUI_programs/Temps/gendata.mat')
mi.nSteps = nSteps;


save('GUI_programs/Temps/gendata.mat','mi','data','-append')


% --- Executes during object creation, after setting all properties.
function enter_nsteps_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
