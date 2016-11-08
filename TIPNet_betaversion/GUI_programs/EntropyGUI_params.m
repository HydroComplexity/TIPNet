function varargout = EntropyGUI_params(varargin)
% ENTROPYGUI_PARAMS MATLAB code for EntropyGUI_params.fig


% Last Modified by GUIDE v2.5 24-May-2016 08:32:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EntropyGUI_params_OpeningFcn, ...
                   'gui_OutputFcn',  @EntropyGUI_params_OutputFcn, ...
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

% --- Executes just before EntropyGUI_params is made visible.
function EntropyGUI_params_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EntropyGUI_params (see VARARGIN)

% Choose default command line output for EntropyGUI_params
handles.output = hObject;

load('GUI_programs/Temps/projectname.mat')
load(fullpathname)

%initialize parameters

set(handles.enter_SigTests,'String',num2str(mi.nTests));
set(handles.enter_nlags,'String',num2str(mi.nlags));

set(handles.parallel_popup,'Value',mi.parallel_opt+1)
set(handles.zero_lag_popup,'Value',mi.ZeroLagOpt+1)

set(handles.self_links_popup,'Value',mi.NoSelfOpt+1);

set(handles.choose_netopt,'String',...
    {'Entropy Only','H(X) and I(X;Y) only','full network'});

set(handles.choose_netopt,'Value',mi.netopt);


%plot lags
axes(handles.plot_lags)
plot(1:mi.nlags,mi.lagvect,'*')
xlabel('Lag number')
ylabel('Lag (timesteps)')


% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = EntropyGUI_params_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;



function enter_SigTests_Callback(hObject, eventdata, handles)

nTests=str2double(get(hObject,'String'));

load('GUI_programs/Temps/projectname.mat')
load(fullpathname)
mi.nTests = nTests;
save(fullpathname,'mi','-append')

% --- Executes during object creation, after setting all properties.
function enter_SigTests_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function enter_nlags_Callback(hObject, eventdata, handles)

nlags=str2double(get(hObject,'String'));
load('GUI_programs/Temps/projectname.mat')
load(fullpathname)
mi.nlags = nlags;
mi.lagvect = 1:nlags;
save(fullpathname,'mi','-append')

%plot lags
axes(handles.plot_lags)
plot(1:mi.nlags,mi.lagvect,'*')
xlabel('Lag number')
ylabel('Lag (timesteps)')

% --- Executes during object creation, after setting all properties.
function enter_nlags_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in button_done.
function button_done_Callback(hObject, eventdata, handles)
% hObject    handle to button_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


close EntropyGUI_params
wind = findobj(0,'Tag','EntropyGUI_main');
set(wind,'Visible','On');


% --- Executes on selection change in choose_netopt.
function choose_netopt_Callback(hObject, eventdata, handles)

load('GUI_programs/Temps/projectname.mat')
load(fullpathname)
mi.netopt = get(hObject,'Value');
save(fullpathname,'mi','-append')

% --- Executes during object creation, after setting all properties.
function choose_netopt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_enter_lagvect.
function button_enter_lagvect_Callback(hObject, eventdata, handles)

addpath('UserData')
filename1 = '/UserData/lagvect.mat';

if exist(filename1,'file') ~= 0
 load(filename1);
 
 load('GUI_programs/Temps/projectname.mat')
 load(fullpathname)
 
 msgbox(sprintf('lag vector file loaded with %d lags between %d and %d timesteps',...
     length(lagvect),min(lagvect),max(lagvect)))
 
 mi.lagvect = lagvect;
 mi.nlags = length(lagvect);
 save(fullpathname,'mi','-append')
 
 set(handles.enter_nlags,'String',num2str(mi.nlags))
 
 %plot lags
axes(handles.plot_lags)
plot(1:mi.nlags,mi.lagvect,'*')
xlabel('Lag number')
ylabel('Lag (timesteps)')
 
else
  msgbox('No file of name "lagvect.mat" found in directory','No file found','error')
end

% --- Executes on selection change in parallel_popup.
function parallel_popup_Callback(hObject, eventdata, handles)

load('GUI_programs/Temps/projectname.mat')
load(fullpathname)
mi.parallel_opt = get(hObject,'Value')-1;
save(fullpathname,'mi','-append')

% --- Executes during object creation, after setting all properties.
function parallel_popup_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in zero_lag_popup.
function zero_lag_popup_Callback(hObject, eventdata, handles)
load('GUI_programs/Temps/projectname.mat')
load(fullpathname)
mi.ZeroLagOpt = get(hObject,'Value')-1;


save(fullpathname,'mi','-append')

% --- Executes during object creation, after setting all properties.
function zero_lag_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zero_lag_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in self_links_popup.
function self_links_popup_Callback(hObject, eventdata, handles)
load('GUI_programs/Temps/projectname.mat')
load(fullpathname)
mi.NoSelfOpt = get(hObject,'Value')-1;
save(fullpathname,'mi','-append')

% --- Executes during object creation, after setting all properties.
function self_links_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to self_links_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
