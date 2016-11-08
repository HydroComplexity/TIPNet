function varargout = EntropyGUI_filter(varargin)
% ENTROPYGUI_FILTER MATLAB code for EntropyGUI_filter.fig

% Last Modified by GUIDE v2.5 11-May-2016 09:36:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EntropyGUI_filter_OpeningFcn, ...
                   'gui_OutputFcn',  @EntropyGUI_filter_OutputFcn, ...
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


% --- Executes just before EntropyGUI_filter is made visible.
function EntropyGUI_filter_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

set(hObject,'toolbar','figure')

load('GUI_programs/Temps/selection.mat');
load('GUI_programs/Temps/projectname.mat')
load(fullpathname) %project file: contains mi structure of model information, data_orig, data

Selected = Selected-1;

dat_orig = mi.DataPrep.dp{Selected}.data_orig;


filter_num = mi.DataPrep.FilterType(Selected+1);
filter_lambda = mi.DataPrep.FilterLambda(Selected+1);

if filter_num ==1
    filtertype = 'hp';
    set(handles.filter_panel,'selectedobject',handles.button_hp)
else
    filtertype = 'lp';
    set(handles.filter_panel,'selectedobject',handles.button_lp)
end



axes(handles.plot_prefilter)
cla(handles.plot_prefilter)
plot(dat_orig)

fy=fft(dat_orig);
py=fy .* conj(fy); % Compute power spectrum
plotrange=2:length(fy)/2;
f=length(dat_orig)./(plotrange-1);

%apply high-pass butterworth filtering
data_filt = ButterFiltFun(filtertype,dat_orig,'lambdac',filter_lambda);
fy=fft(data_filt);
py2=fy .* conj(fy); % Compute new power spectrum

axes(handles.plot_powerspectra)
cla(handles.plot_powerspectra)
semilogy(1./f,real(py(plotrange)),'-.b','LineWidth',2)
hold on
%plot post-filtering power spectra
semilogy(1./f,real(py2(plotrange)),'-.r','LineWidth',2)
legend('original','filtered','Location','EastOutside')
xlabel('frequency')
ylabel('power')

axes(handles.plot_postfilter)
cla(handles.plot_postfilter)
plot(data_filt)
ylim([min(data_filt) max(data_filt)])


mi.DataPrep.dp{Selected}.data_process2 = data_filt;
mi.DataPrep.dp{Selected}.data_process1 = data_filt;

%update selection in "all nodes" matrix
mi.DataPrep.alldata_process1(:,Selected)=data_filt;
mi.DataPrep.alldata_process2(:,Selected)=data_filt;

handles.mi = mi;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = EntropyGUI_filter_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;



% --- Executes on button press in button_donefilt.
function button_donefilt_Callback(hObject, eventdata, handles)

load('GUI_programs/Temps/selection.mat');
load('GUI_programs/Temps/projectname.mat')
%load(fullpathname) %project file: contains mi structure of model information, data_orig, data

mi = handles.mi;
save(fullpathname,'mi','-append');

close(EntropyGUI_filter);
activefig = findobj(0,'Tag','EntropyGUI_prep');
set(activefig,'Visible','On');


% --- Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)
%reset processed data to original data
load('GUI_programs/Temps/projectname.mat')
load('GUI_programs/Temps/selection.mat');
Selected = Selected -1;

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

close EntropyGUI_filter


% --- Executes when selected object is changed in filter_panel.
function filter_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in filter_panel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

mi = handles.mi;
load('GUI_programs/Temps/selection.mat');
Selected = Selected -1;

Choice=get(eventdata.NewValue,'Tag'); % Get Tag of selected object

switch Choice % Get Tag of selected object.
    case 'button_hp'  
        mi.DataPrep.FilterType(Selected+1) = 1;
    case 'button_lp'
        mi.DataPrep.FilterType(Selected+1) = 2;
end

handles.mi = mi;

guidata(hObject,handles)


% --- Executes on button press in button_plot.
function button_plot_Callback(hObject, eventdata, handles)

mi = handles.mi;

load('GUI_programs/Temps/projectname.mat')
load('GUI_programs/Temps/selection.mat');
Selected = Selected -1;

dat_orig = mi.DataPrep.dp{Selected}.data_orig;

filter_num = mi.DataPrep.FilterType(Selected+1);
filter_lambda = mi.DataPrep.FilterLambda(Selected+1);

if filter_num ==1
    filtertype = 'hp';
else
    filtertype = 'lp';
end

axes(handles.plot_prefilter)
cla(handles.plot_prefilter)
plot(dat_orig)

fy=fft(dat_orig);
py=fy .* conj(fy); % Compute power spectrum
plotrange=2:length(fy)/2;
f=length(dat_orig)./(plotrange-1);

%apply high-pass butterworth filtering
data_filt = ButterFiltFun(filtertype,dat_orig,'lambdac',filter_lambda);
fy=fft(data_filt);
py2=fy .* conj(fy); % Compute new power spectrum

axes(handles.plot_powerspectra)
cla(handles.plot_powerspectra)
semilogy(1./f,real(py(plotrange)),'-.b','LineWidth',2)
hold on
%plot post-filtering power spectra
semilogy(1./f,real(py2(plotrange)),'-.r','LineWidth',2)
legend('original','filtered','Location','EastOutside')
xlabel('frequency')
ylabel('power')

axes(handles.plot_postfilter)
cla(handles.plot_postfilter)
plot(data_filt)
ylim([min(data_filt) max(data_filt)])

mi.DataPrep.dp{Selected}.data_process2 = data_filt;
mi.DataPrep.dp{Selected}.data_process1 = data_filt;

%update selection in "all nodes" matrix
mi.DataPrep.alldata_process1(:,Selected)=data_filt;
mi.DataPrep.alldata_process2(:,Selected)=data_filt;

handles.mi = mi;
save(fullpathname,'mi','-append')

% Update handles structure
guidata(hObject, handles);
