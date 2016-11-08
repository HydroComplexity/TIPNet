function varargout = EntropyGUI_preprocess(varargin)
% ENTROPYGUI_PREPROCESS MATLAB code for EntropyGUI_preprocess.fig
% Last Modified by GUIDE v2.5 10-May-2016 21:09:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EntropyGUI_preprocess_OpeningFcn, ...
                   'gui_OutputFcn',  @EntropyGUI_preprocess_OutputFcn, ...
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

% --- Executes just before EntropyGUI_preprocess is made visible.
function EntropyGUI_preprocess_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

set(hObject,'toolbar','figure')

%populate listbox to have choices of all individual nodes or all nodes
%together

load('GUI_programs/Temps/projectname.mat')
load(fullpathname) %project file: contains mi structure of model information, data_orig, data

handles.mi = mi;
handles.netdata = netdata;
handles.plotlim = min(10000,mi.nSteps);
handles.varnames = varnames;
plotlim = handles.plotlim;

choices = ['All Nodes',varnames];
set(handles.select_nodes,'String',choices)

Selected=1;
set(handles.select_nodes,'Value',1)
save('GUI_programs/Temps/selection.mat','Selected');

 %can only filter one node at a time
set(handles.button_filter,'Enable','Off')
set(handles.button_anomaly,'Enable','Off')

set(handles.enter_segsteps,'String',num2str(mi.segsteps));

if mi.DataPrep.FilterOpt(1)==1
set(handles.select_filtering,'selectedobject',handles.button_none)
elseif mi.DataPrep.FilterOpt(1)==2
set(handles.select_filtering,'selectedobject',handles.button_anomaly)
elseif mi.DataPrep.FilterOpt(1)==3
set(handles.select_filtering,'selectedobject',handles.button_filter)
elseif mi.DataPrep.FilterOpt(1)==4
set(handles.select_filtering,'selectedobject',handles.button_increment)
else
set(handles.select_filtering,'selectedobject',handles.button_log)    
end

if mi.DataPrep.NormOpt(1)==1
set(handles.check_normalize,'Value',1)
end

if mi.DataPrep.OutlierOpt(1)==1
set(handles.check_removeoutliers,'Value',1)
end

%%%%%%%%%%%%%%%%%%%%%%%% Initialize Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.plot_processing)
cla(handles.plot_processing,'reset')
hold on
nvars = mi.nvars;
cvect = jet(mi.nvars);
for n =1:nvars
plot(mi.DataPrep.dp{n}.data_process2(1:plotlim),'Color',cvect(n,:))
end
xlim([1 plotlim])


axes(handles.plot_orig)
cla(handles.plot_orig,'reset')
hold on
for n =1:nvars
plot(mi.DataPrep.dp{n}.data_orig(1:plotlim),'Color',cvect(n,:))
end
xlim([1 plotlim])

axes(handles.LegSpace)
hold on
set(gca,'Visible','off')
for n =1:nvars
plot(mi.DataPrep.dp{n}.data_orig(1:plotlim),'Color',cvect(n,:))
end
xlim([1 2])
leg=legend(varnames,'Location','West');
pos = get(leg,'Position');
set(leg,'Position',pos +[-5 0 0 0])


axes(handles.plot_all)
cla(handles.plot_all,'reset')
hold on
for n =1:nvars
plot(mi.DataPrep.dp{n}.data_process2(1:plotlim),'Color',cvect(n,:))
end
xlim([1 plotlim])


x =0;
minval = min(min(mi.DataPrep.alldata_process2));
maxval = max(max(mi.DataPrep.alldata_process2));
for i =1:mi.nsegs
    x = x+mi.segsteps;
    line([x x],[minval maxval],'LineWidth',2,'Color','k')
end
xlim([1 plotlim])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%disable anomaly option if too few data points
if size(mi.DataPrep.dp{1}.data_orig,1)<200
set(handles.button_anomaly,'Enable','off')
end

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = EntropyGUI_preprocess_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


% --- Executes on selection change in select_nodes.
function select_nodes_Callback(hObject, eventdata, handles)

Selected = get(hObject,'Value'); %first index is "all nodes", then nodes in order

mi = handles.mi;

check_norm = mi.DataPrep.NormOpt(Selected);
check_out = mi.DataPrep.OutlierOpt(Selected);

set(handles.check_removeoutliers,'Value',check_out);
set(handles.check_normalize,'Value',check_norm);

%set(handles.select_filtering,'selectedobject',handles.button_none)
if mi.DataPrep.FilterOpt(Selected)==1
    set(handles.select_filtering,'selectedobject',handles.button_none)
elseif mi.DataPrep.FilterOpt(Selected)==2
    set(handles.select_filtering,'selectedobject',handles.button_anomaly)
elseif mi.DataPrep.FilterOpt(Selected)==3
     set(handles.select_filtering,'selectedobject',handles.button_filter)
elseif mi.DataPrep.FilterOpt(Selected)==4
    set(handles.select_filtering,'selectedobject',handles.button_increment)
else
    set(handles.select_filtering,'selectedobject',handles.button_log)
end

if Selected==1 %can only filter one node at a time
    set(handles.button_filter,'Enable','Off')
    set(handles.button_anomaly,'Enable','Off')
else
    set(handles.button_filter,'Enable','on')
    set(handles.button_anomaly,'Enable','on')
end

save('GUI_programs/Temps/selection.mat','Selected')

set(handles.button_plot,'backgroundcolor','y')

% --- Executes during object creation, after setting all properties.
function select_nodes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in reset_preprocessing.
function reset_preprocessing_Callback(hObject, eventdata, handles)
% reset any altered data to orginally loaded data for all nodes

%load('GUI_programs/Temps/selection.mat');
load('GUI_programs/Temps/projectname.mat')
%load(fullpathname) %project file: contains mi structure of model information, data_orig, data

mi = handles.mi;
varnames = handles.varnames;
plotlim = handles.plotlim;

    for i = 1:mi.nvars
    mi.DataPrep.dp{i}.data_process2=mi.DataPrep.dp{i}.data_orig;
    mi.DataPrep.dp{i}.data_process1=mi.DataPrep.dp{i}.data_orig;
    end
    
    mi.DataPrep.alldata_process1 = mi.DataPrep.alldata_orig;
    mi.DataPrep.alldata_process2 = mi.DataPrep.alldata_orig;
    
    mi.DataPrep.FilterOpt(:)=1;
    mi.DataPrep.OutlierOpt(:)=0;
    mi.DataPrep.NormOpt(:)=0;

mi.nsegs=1;
mi.segsteps = mi.nSteps;

%%%%%%%%%%%%%%%%%%%%%%%% Reset Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.plot_processing)
cla(handles.plot_processing,'reset')
hold on
cvect = jet(mi.nvars);
nvars = mi.nvars;
for n =1:nvars
plot(mi.DataPrep.dp{n}.data_process2(1:plotlim),'Color',cvect(n,:))
end
xlim([1 plotlim])


axes(handles.plot_orig)
cla(handles.plot_orig,'reset')
hold on
for n =1:nvars
plot(mi.DataPrep.dp{n}.data_orig(1:plotlim),'Color',cvect(n,:))
end
xlim([1 plotlim])


axes(handles.plot_all)
cla(handles.plot_all,'reset')
hold on
for n =1:nvars
plot(mi.DataPrep.dp{n}.data_process2(1:plotlim),'Color',cvect(n,:))
end
xlim([1 plotlim])


x =0;
minval = min(min(mi.DataPrep.alldata_process2));
maxval = max(max(mi.DataPrep.alldata_process2));
for i =1:mi.nsegs
    x = x+mi.segsteps;
    line([x x],[minval maxval],'LineWidth',2,'Color','k')
end
xlim([1 plotlim])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set node selection to "all nodes"
set(handles.select_nodes,'Value',1)
Selected = 1;
save('GUI_programs/Temps/selection.mat','Selected');

set(handles.select_filtering,'selectedobject',handles.button_none)
set(handles.check_removeoutliers,'Value',0);
set(handles.check_normalize,'Value',0);
set(handles.enter_segsteps,'String',num2str(mi.segsteps));

%revert to single segment of data
netdata=[];
netdata{1}.data_orig = mi.DataPrep.alldata_orig;
netdata{1}.data_process1 = mi.DataPrep.alldata_orig;
netdata{1}.data_process2 = mi.DataPrep.alldata_orig;

handles.mi = mi;
handles.netdata = netdata;

save(fullpathname,'mi','netdata','-append')
guidata(hObject,handles)
%save(fullpathname,'mi','-append');

% --- Executes on button press in check_removeoutliers.
function check_removeoutliers_Callback(hObject, eventdata, handles)
% remove outliers (set to upper and lower bounds)
check_rmout = get(handles.check_removeoutliers,'Value');
check_norm = get(handles.check_normalize,'Value');


mi = handles.mi;

load('GUI_programs/Temps/projectname.mat')
load('GUI_programs/Temps/selection.mat');

mi.DataPrep.OutlierOpt(Selected)=check_rmout;

if Selected==1 && mi.DataPrep.OutlierOpt(1)==1
    mi.DataPrep.OutlierOpt(:)=1;
elseif Selected ==1 && mi.DataPrep.OutlierOpt(1)==0
    mi.DataPrep.OutlierOpt(:)=0;
end

%update
if Selected>1 %if an individual node is changed
    mi.DataPrep.dp{Selected-1}.data_process2 = ProcessFunction(mi.DataPrep.dp{Selected-1}.data_process1,check_rmout,check_norm);
    mi.DataPrep.alldata_process2(:,Selected-1)=mi.DataPrep.dp{Selected-1}.data_process2;
elseif Selected ==1 %all nodes changed
    
    mi.DataPrep.alldata_process2 = ProcessFunction(mi.DataPrep.alldata_process1,check_rmout,check_norm);
    for i = 1:mi.nvars    
    mi.DataPrep.dp{i}.data_process2=mi.DataPrep.alldata_process2(:,i);
    end
end


handles.mi = mi;
guidata(hObject,handles);
save(fullpathname,'mi','-append')

set(handles.button_plot,'backgroundcolor','y')

% --- Executes on button press in check_normalize.
function check_normalize_Callback(hObject, eventdata, handles)
% normalize values between 0 and 1
check_rmout = get(handles.check_removeoutliers,'Value');
check_norm = get(handles.check_normalize,'Value');

load('GUI_programs/Temps/projectname.mat')
load('GUI_programs/Temps/selection.mat');

mi = handles.mi;

mi.DataPrep.NormOpt(Selected)=check_norm;

if Selected==1 && mi.DataPrep.NormOpt(1)==1
    mi.DataPrep.NormOpt(:)=1;
elseif Selected ==1 && mi.DataPrep.NormOpt(1)==0
    mi.DataPrep.NormOpt(:)=0;
end

%update
if Selected>1 %if an individual node is changed
    mi.DataPrep.dp{Selected-1}.data_process2 = ProcessFunction(mi.DataPrep.dp{Selected-1}.data_process1,check_rmout,check_norm);
    mi.DataPrep.alldata_process2(:,Selected-1)=mi.DataPrep.dp{Selected-1}.data_process2;
elseif Selected ==1 %all nodes changed
    
    mi.DataPrep.alldata_process2 = ProcessFunction(mi.DataPrep.alldata_process1,check_rmout,check_norm);
    for i = 1:mi.nvars    
    mi.DataPrep.dp{i}.data_process2=mi.DataPrep.alldata_process2(:,i);
    end
end

save(fullpathname,'mi','-append')
set(handles.button_plot,'backgroundcolor','y')

handles.mi = mi;
guidata(hObject,handles)

function enter_segsteps_Callback(hObject, eventdata, handles)

%load('GUI_programs/Temps/projectname.mat')
%load(fullpathname) %project file: contains mi structure of model information, data_orig, data
load('GUI_programs/Temps/projectname.mat')

mi = handles.mi;

mi.segsteps = str2double(get(hObject,'String'));
mi.nsegs = floor(mi.nSteps./mi.segsteps);

if mi.segsteps>mi.nSteps
    msgbox('Segment size must be lower than total data points')
    mi.segsteps = mi.nSteps;
else
    %save(fullpathname,'mi','-append')
    set(handles.button_plot,'backgroundcolor','y')
end

save(fullpathname,'mi','-append')

handles.mi = mi;
guidata(hObject,handles)
  

function enter_segsteps_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in select_filtering.
function select_filtering_SelectionChangeFcn(hObject, eventdata, handles)

load('GUI_programs/Temps/projectname.mat')
load('GUI_programs/Temps/selection.mat');
mi = handles.mi;


Choice=get(eventdata.NewValue,'Tag'); % Get Tag of selected object
check_rmout = get(handles.check_removeoutliers,'Value');
check_norm = get(handles.check_normalize,'Value');

%save(fullpathname,'mi','-append') %save any other changes made to mi
load(fullpathname) %project file: contains mi structure of model information, data_orig, data

%reset selection or all nodes to original value
if Selected >1
mi.DataPrep.dp{Selected-1}.data_process1=mi.DataPrep.dp{Selected-1}.data_orig;
mi.DataPrep.dp{Selected-1}.data_process2=mi.DataPrep.dp{Selected-1}.data_orig;
else
mi.DataPrep.alldata_process2 = mi.DataPrep.alldata_orig;
end

save(fullpathname,'mi','-append')

switch Choice % Get Tag of selected object.
    case 'button_none'      %no action, revert to original data
        mi.DataPrep.FilterOpt(Selected)=1;            
    case 'button_filter'    %butterworth filtering, enter new fig window
        mi.DataPrep.FilterOpt(Selected)=3;
        save(fullpathname,'mi','-append')
        activefig = findobj(0,'Tag','EntropyGUI_prep');
        set(activefig,'Visible','Off');
        EntropyGUI_filter; 
        waitfor(activefig,'Visible')
        load(fullpathname)
    case 'button_anomaly'   %take anomaly to filter cycle, enter new fig window
   
    mi.DataPrep.FilterOpt(Selected)=2;
    save(fullpathname,'mi','-append')
    activefig = findobj(0,'Tag','EntropyGUI_prep');
    set(activefig,'Visible','Off');
    EntropyGUI_prep_anomchoose;
    waitfor(activefig,'Visible')
    load(fullpathname)
    
    case 'button_log'       %take log10 of data
        if Selected >1
            dat = log10(mi.DataPrep.dp{Selected-1}.data_orig);
            flag = zeros(size(dat));
            flag(dat<-5)=1;
            mindat = min(dat(flag==0));
            dat(dat<mindat)=mindat;
            mi.DataPrep.dp{Selected-1}.data_process1 = dat;
            mi.DataPrep.dp{Selected-1}.data_process2 = dat;
        else
            dat = log10(mi.DataPrep.alldata_orig);
            flag = zeros(size(dat));
            flag(dat<-5)=1;
            mindat = min(dat(flag==0));
            dat(dat<mindat)=mindat;
            mi.DataPrep.alldata_process1 = dat;
            mi.DataPrep.alldata_process2 = dat;
        end
        mi.DataPrep.FilterOpt(Selected)=5;
    case 'button_increment'  %take increment of data
        
        if Selected >1
            X = mi.DataPrep.dp{Selected-1}.data_orig(2:end,:)-mi.DataPrep.dp{Selected-1}.data_orig(1:end-1,:);
            mi.DataPrep.dp{Selected-1}.data_process1 = [0; X];
            mi.DataPrep.dp{Selected-1}.data_process2 = [0; X];
        else
            
            X = mi.DataPrep.alldata_orig(2:end,:)-mi.DataPrep.alldata_orig(1:end-1,:);
            mi.DataPrep.alldata_process2 = [zeros(size(mi.DataPrep.alldata_orig,2)); X];
            mi.DataPrep.alldata_process1 = [zeros(size(mi.DataPrep.alldata_orig,2)); X];
        end
     mi.DataPrep.FilterOpt(Selected)=4;
end

if Selected==1
 mi.DataPrep.FilterOpt(:)=mi.DataPrep.FilterOpt(1);
end

%update data_process1 and data_process2
if Selected>1 %if an individual node is reset
    
    mi.DataPrep.dp{Selected-1}.data_process2 = ProcessFunction(mi.DataPrep.dp{Selected-1}.data_process1,check_rmout,check_norm);
    mi.DataPrep.alldata_process2(:,Selected-1) = mi.DataPrep.dp{Selected-1}.data_process2;

elseif Selected ==1
    
    mi.DataPrep.alldata_process2 = ProcessFunction(mi.DataPrep.alldata_process1,check_rmout,check_norm);
    for i = 1:mi.nvars
    mi.DataPrep.dp{i}.data_process2=mi.DataPrep.alldata_process2(:,i);
    end
end

handles.mi = mi;
set(handles.button_plot,'backgroundcolor','y')
save(fullpathname,'mi','-append');

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function select_filtering_CreateFcn(hObject, eventdata, handles)

% --- Executes on button press in button_done.
function button_done_Callback(hObject, eventdata, handles)

load('GUI_programs/Temps/projectname.mat')
load(fullpathname) %project file: contains mi structure of model information, data_orig, data

if mi.segsteps <= mi.nSteps
  
  netdata=[];
  
  mi.nsegs = floor(mi.nSteps/mi.segsteps);

  start = 1;
  stop = mi.segsteps;
  for i = 1:mi.nsegs
      netdata{i}.data_orig = mi.DataPrep.alldata_orig(start:stop,:);
      netdata{i}.data_process1 = mi.DataPrep.alldata_process1(start:stop,:);
      netdata{i}.data_process2 = mi.DataPrep.alldata_process2(start:stop,:);
      start = stop+1;
      stop = start+mi.segsteps-1;
  end
   
end

%reset range if data has been altered
mi.Range = [min(mi.DataPrep.alldata_process2); max(mi.DataPrep.alldata_process2)];


%set h (pdf smoothing parameter) if data has been altered
for i=1:mi.nvars
    
  ndata = sum(mi.DataPrep.dp{i}.data_process2~=0); %number non-zero values
  nsteps = mi.segsteps .*(ndata./mi.nSteps); %timesteps per segment * percentage of non-zeros
  
  if nsteps>=1
  mi.KDEparams.h1_seg_orig(i)= 1.06 * nsteps^(-1/5) .* var(mi.DataPrep.dp{i}.data_process2); 
  mi.KDEparams.h2_seg_orig(i)= 1.77 * nsteps^(-1/6) .* var(mi.DataPrep.dp{i}.data_process2); 
  mi.KDEparams.h3_seg_orig(i) = 2.78 * nsteps^(-1/7) .* var(mi.DataPrep.dp{i}.data_process2); 
  else                          %all zero values, no pdf anyway
  mi.KDEparams.h1_seg_orig(i) = 0; 
  mi.KDEparams.h2_seg_orig(i) = 0; 
  mi.KDEparams.h3_seg_orig(i) = 0; 
  end
  
end

mi.KDEparams.h1 = mi.KDEparams.h1_seg_orig;
mi.KDEparams.h2 = mi.KDEparams.h2_seg_orig;
mi.KDEparams.h3 = mi.KDEparams.h3_seg_orig;

save(fullpathname,'mi','netdata','-append')

close EntropyGUI_preprocess

wind = findobj(0,'Tag','EntropyGUI_main');
set(wind,'Visible','On')


% --- Executes on button press in button_plot.
function button_plot_Callback(hObject, eventdata, handles)

set(handles.button_plot,'backgroundcolor',[0.9412    0.9412    0.9412])

mi = handles.mi;
varnames = handles.varnames;
cvect1 = jet(mi.nvars);


load('GUI_programs/Temps/selection.mat');
cvect2 = ones(size(cvect1));
if Selected==1
  cvect2 = cvect1;
else
    cvect2(Selected-1,:)=cvect1(Selected-1,:);
end

plotlim = handles.plotlim;

nvars =mi.nvars;

%plots
axes(handles.plot_processing)
cla(handles.plot_processing,'reset')
hold on
cvect = jet(mi.nvars);

for n =1:nvars
plot(mi.DataPrep.dp{n}.data_process2(1:plotlim),'Color',cvect2(n,:))
end

if Selected>1
plot(mi.DataPrep.dp{Selected-1}.data_process2(1:plotlim),'Color',cvect2(Selected-1,:))
ylim([min(mi.DataPrep.dp{Selected-1}.data_process2)-.001 max(mi.DataPrep.dp{Selected-1}.data_process2)+.001])
end

xlim([1 plotlim])


axes(handles.plot_orig)
cla(handles.plot_orig,'reset')
hold on
for n =1:nvars
plot(mi.DataPrep.dp{n}.data_orig(1:plotlim),'Color',cvect2(n,:))
end
if Selected>1
plot(mi.DataPrep.dp{Selected-1}.data_orig(1:plotlim),'Color',cvect2(Selected-1,:))
ylim([min(mi.DataPrep.dp{Selected-1}.data_orig)-.001 max(mi.DataPrep.dp{Selected-1}.data_orig)+.001])
end
xlim([1 plotlim])


axes(handles.plot_all)
cla(handles.plot_all,'reset')
hold on
for n =1:nvars
plot(mi.DataPrep.dp{n}.data_process2(1:plotlim),'Color',cvect1(n,:))
end
xlim([1 plotlim])



x =0;
minval = min(min(mi.DataPrep.alldata_process2));
maxval = max(max(mi.DataPrep.alldata_process2));
for i =1:mi.nsegs
    x = x+mi.segsteps;
    line([x x],[minval maxval],'LineWidth',2,'Color','k')
end
xlim([1 plotlim])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
