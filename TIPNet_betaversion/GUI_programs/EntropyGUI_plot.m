function varargout = EntropyGUI_plot(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EntropyGUI_plot_OpeningFcn, ...
                   'gui_OutputFcn',  @EntropyGUI_plot_OutputFcn, ...
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


% --- Executes just before EntropyGUI_plot is made visible.
function EntropyGUI_plot_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for EntropyGUI_plot
handles.output = hObject;


set(hObject,'toolbar','figure')

load('GUI_programs/Temps/projectname.mat')
load(fullpathname) %project file: contains mi structure of model information, data_orig, data

%populate popup menu_segment to choose between "averages" or individual
%segments
if mi.nsegs>1
choices = 'averages';
for i=1:mi.nsegs
   segs{i}=i; 
end
choices =[choices,segs];
set(handles.menu_segment,'String',choices)
else %only one segment: don't show choice or window timeseries plot
set(handles.menu_segment,'Visible','Off') 
set(handles.text_chooseseg,'Visible','Off')
end

set(handles.button_plot,'backgroundcolor',[0.9412    0.9412    0.9412])

%populate popup menu_lag to choose between max I value for all lags or
%individual lag times
lagvect=mi.lagvect;
nlags = mi.nlags;

if mi.ZeroLagOpt==1
   lagvect =[0 lagvect]; 
   nlags=nlags+1;
end

for i =1:nlags
   lags{i}=lagvect(i); 
end

choices = ['dominant links','zero lag',lags];
set(handles.menu_lags,'String',choices);

%populate popup menu_transnode and menu_recnode to choose nodes or all
%links
choices = ['all nodes', cellstr(varnames)];

set(handles.menu_transnode,'String',choices)
set(handles.menu_recnode,'String',choices)

%populate popup menu_plottype for types of network measures to plot
choices = {'Mutual Information I','Transfer Entropy T_E','T/I',...
    'Unique Info','Link Entropy H_L','Synergy S','Redundancy R'};
set(handles.menu_plottype,'String',choices)

%Initialize Plots for lagged MI (dominant lags, all segments averaged)
if mi.nsegs>1
I_mut = AllStats.I_normbyH;
I_lag = AllStats.I_lag;
else
I_mut = entropy{1}.I_dom;
I_lag = entropy{1}.I_dom_tau;
end

titletext = '';
uplim = .0001;
bothways = 1;

Links = I_mut;
Weights = I_lag;
A = 1;
WeightRange = [0 max(lagvect) length(lagvect)];

for i =1:mi.nvars
    nodesize(i) = max(.01,Links(i,i));
end

colorlabel = 'lag time';
handles.pos = get(handles.plot_circle,'position');
cla(handles.plot_circle);
fh=handles.plot_circle;

PlotFunGUI(fh, Links,Weights, WeightRange, nodesize, varnames,titletext,A,uplim,bothways);

axes(handles.plot_legend);
fh=handles.plot_legend;
PlotLegendGUI(fh,lagvect,colorlabel)

cla(handles.plot_windows,'reset');
axes(handles.plot_windows)
hold on;

for i =1:mi.nsegs
    SegI(i) = mean(mean(entropy{i}.I_dom));
    SegTI(i) = mean(mean(entropy{i}.TI_T));
    SegTE(i) = mean(mean(entropy{i}.TE_T_normbyItot));
    SegIinst(i) = mean(mean(entropy{i}.I_normbyH(1,:,:)));
    SegS(i) = mean(mean(entropy{i}.S_T));
    SegR(i) = mean(mean(entropy{i}.R_T));
    SegU(i) = mean(mean(entropy{i}.U_T));
end

maxval = max(max([SegI SegTI SegTE SegIinst SegS SegR SegU]));

Lwid = ones(1,7);
Lwid(1)=3;
Msize = ones(1,7)+3;
Msize(1)=6;

plot(1:mi.nsegs,SegI,'-o','Color','b','LineWidth',Lwid(1),'MarkerSize',Msize(1))
plot(1:mi.nsegs,SegIinst,'-o','Color','g','LineWidth',Lwid(1),'MarkerSize',Msize(1))
plot(1:mi.nsegs,SegTE,'-o','Color','m','LineWidth',Lwid(2),'MarkerSize',Msize(2))
plot(1:mi.nsegs,SegTI,'-o','Color','k','LineWidth',Lwid(3),'MarkerSize',Msize(3))
%plot(1:mi.nsegs,SegU,'-o','Color',[0 .6 0],'LineWidth',Lwid(4),'MarkerSize',Msize(4))
%plot(1:mi.nsegs,SegS,'-o','Color',[.8 0 0],'LineWidth',Lwid(6),'MarkerSize',Msize(6))
%plot(1:mi.nsegs,SegR,'-o','Color','c','LineWidth',Lwid(7),'MarkerSize',Msize(7))


legend('Lagged I','T/I','T_E','Inst I','Location','East')
%xlabel('Segment')
title('Average Quantities (normalized)')
ylim([0 maxval])
xlim([1 max(mi.nsegs,1.01)])
ylabel('bits/bit')

%bar plot
cla(handles.plot_SRU,'reset');
axes(handles.plot_SRU)
hold on;

mc_U = [78 45 109] ./ 255;
mc_S = [238 240 2] ./ 255;
mc_R = [.7 0 0];

myC = [mc_S; mc_R; mc_U];

y = [SegS; SegR; SegU]';
H = bar(1:mi.nsegs, y,'stacked');
for k=1:3
  set(H(k),'facecolor',myC(k,:))
end

xlim([1 max(mi.nsegs,1.01)])
xlabel('Segment')
ylabel('bits')
legend('S','R','U','Location','East')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EntropyGUI_plot wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = EntropyGUI_plot_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on selection change in menu_segment.
function menu_segment_Callback(hObject, eventdata, handles)

%if certain segment selected while menu_plottype value = 5 (HLinks)
%error message and set measure to 1 (information)
type = get(handles.menu_plottype,'Value');
seg = get(handles.menu_segment,'Value');
if seg>1 && type==5
    msgbox('Warning: no H(Link) for single segment, reverting to mutual information',...
        'error','error')
    set(handles.menu_plottype,'Value',1)
end

set(handles.button_plot,'backgroundcolor','y')

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function menu_segment_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in menu_lags.
function menu_lags_Callback(hObject, eventdata, handles)

%if lag>1 (part of lagvect or zero lag), exclude T/I,U, and H(Link) as
%types

set(handles.button_plot,'backgroundcolor','y')

lag = get(handles.menu_lags,'Value');
type = get(handles.menu_plottype,'Value');
if lag>1 && type > 2 %T/I =3, U=4, H=5, S=6, R=7
    msgbox('Warning: T/I, U, H(Link), R, S only for dominant time lags, reverting to mutual information',...
        'error','error')
    set(handles.menu_plottype,'Value',1)   
end

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function menu_lags_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in button_plot.
function button_plot_Callback(hObject, eventdata, handles)

set(handles.button_plot,'backgroundcolor',[0.9412    0.9412    0.9412])

load('GUI_programs/Temps/projectname.mat')
load(fullpathname) %project file: contains mi structure of model information, data_orig, data

nvars = mi.nvars;
lagvect = mi.lagvect;

if mi.ZeroLagOpt ==1
    lagvect =[0 lagvect];
end

time_pick = get(handles.menu_lags,'Value');

seg = get(handles.menu_segment,'Value');
transnode = get(handles.menu_transnode,'Value');
recnode = get(handles.menu_recnode,'Value');
type = get(handles.menu_plottype,'Value');

%make line width thicker for chosen metric for multiple segments
Lwid = ones(1,7);
Lwid(type)=3;
Msize = ones(1,7)+3;
Msize(type)=6;

if transnode>1 && recnode > 1 %only show segment time series for pair of nodes
     
%find nodepair index trans and rec
Trans = transnode-1;
Rec = recnode-1;

%plot timeseries for selected node pair
    cla(handles.plot_windows,'reset');
    axes(handles.plot_windows)
    hold on;
    for i =1:mi.nsegs
        SegI(i) = entropy{i}.I_dom(Trans,Rec);
        SegTI(i) = entropy{i}.TI_T(Trans,Rec);
        SegTE(i) = entropy{i}.TE_T_normbyItot(Trans,Rec);
        SegIinst(i) = entropy{i}.I_normbyH(1,Trans,Rec);
        SegS(i) = entropy{i}.S_T(Trans,Rec);
        SegR(i) = entropy{i}.R_T(Trans,Rec);
        SegU(i) = entropy{i}.U_T(Trans,Rec);
    end
    
elseif transnode==1 && recnode>1  %plot recieved information to one node from all
    %find nodepair index trans and rec
    Rec = recnode-1;

    %plot timeseries for selected node pair
    cla(handles.plot_windows,'reset');
    axes(handles.plot_windows)
    hold on;
    for i =1:mi.nsegs
        SegI(i) = mean(entropy{i}.I_dom(:,Rec));
        SegTI(i) = mean(entropy{i}.TI_T(:,Rec));
        SegTE(i) = mean(entropy{i}.TE_T_normbyItot(:,Rec));
        SegIinst(i) = mean(entropy{i}.I_normbyH(1,:,Rec));
        SegS(i) = mean(entropy{i}.S_T(:,Rec));
        SegR(i) = mean(entropy{i}.R_T(:,Rec));
        SegU(i) = mean(entropy{i}.U_T(:,Rec));
    end
    

elseif recnode == 1 && transnode>1 %plot transmitted info from one node to all
        %find nodepair index trans and rec
    Trans = transnode-1;

    %plot timeseries for selected node pair
    cla(handles.plot_windows,'reset');
    axes(handles.plot_windows)
    hold on;
    for i =1:mi.nsegs
        SegI(i) = mean(entropy{i}.I_dom(Trans,:));
        SegTI(i) = mean(entropy{i}.TI_T(Trans,:));
        SegTE(i) = mean(entropy{i}.TE_T_normbyItot(Trans,:));
        SegIinst(i) = mean(entropy{i}.I_normbyH(1,Trans,:));
        SegS(i) = mean(entropy{i}.S_T(Trans,:));
        SegR(i) = mean(entropy{i}.R_T(Trans,:));
        SegU(i) = mean(entropy{i}.U_T(Trans,:));
    end
    

    
else %plot averages for all values
    
    cla(handles.plot_windows,'reset');
    axes(handles.plot_windows)
    hold on;
    for i =1:mi.nsegs
        SegI(i) = mean(mean(entropy{i}.I_dom));
        SegTI(i) = mean(mean(entropy{i}.TI_T));
        SegTE(i) = mean(mean(entropy{i}.TE_T_normbyItot));
        SegIinst(i) = mean(mean(entropy{i}.I_normbyH(1,:,:)));
        SegS(i) = mean(mean(entropy{i}.S_T));
        SegR(i) = mean(mean(entropy{i}.R_T));
        SegU(i) = mean(mean(entropy{i}.U_T));
    end
     
end

%plot
plot(1:mi.nsegs,SegI,'-o','Color','b','LineWidth',Lwid(1),'MarkerSize',Msize(1))
plot(1:mi.nsegs,SegIinst,'-o','Color','g','LineWidth',Lwid(1),'MarkerSize',Msize(1))
plot(1:mi.nsegs,SegTE,'-o','Color','m','LineWidth',Lwid(2),'MarkerSize',Msize(2))
plot(1:mi.nsegs,SegTI,'-o','Color','k','LineWidth',Lwid(3),'MarkerSize',Msize(3))
%plot(1:mi.nsegs,SegU,'-o','Color',[0 .6 0],'LineWidth',Lwid(4),'MarkerSize',Msize(4))
%plot(1:mi.nsegs,SegS,'-o','Color',[.8 0 0],'LineWidth',Lwid(6),'MarkerSize',Msize(6))
%plot(1:mi.nsegs,SegR,'-o','Color','c','LineWidth',Lwid(7),'MarkerSize',Msize(7))


title('Average Quantities')

maxval = max(max([SegI SegTI SegTE SegIinst .001]));

ylim([0 maxval])
xlim([1 max(mi.nsegs,1.01)])

ylabel('bits/bit')
legend('Lagged I','Inst I','T_E','T/I','Location','East')

%bar plot
cla(handles.plot_SRU,'reset');
axes(handles.plot_SRU)
hold on;

mc_U = [78 45 109] ./ 255;
mc_S = [238 240 2] ./ 255;
mc_R = [.7 0 0];

myC = [mc_S; mc_R; mc_U];

y = [SegS; SegR; SegU]';
H = bar(1:mi.nsegs, y,'stacked');
for k=1:3
  set(H(k),'facecolor',myC(k,:))
end

xlim([1 max(mi.nsegs,1.01)])
xlabel('Segment')
ylabel('bits')
legend('S','R','U','Location','East')


if mi.nsegs==1
    AllStats = entropy{1};
    AllStats.I_normbyH = AllStats.I_dom;
    AllStats.I_lag = AllStats.I_dom_tau;
    AllStats.TI = AllStats.TI_T;
    AllStats.TE = AllStats.TE_T_normbyItot;
    AllStats.U = AllStats.U_T_normbyItot;
    AllStats.R = AllStats.R_T_normbyItot;
    AllStats.S = AllStats.S_T_normbyItot;
    seg=1;
end

TE = zeros(nvars); TI = zeros(nvars); U = zeros(nvars); H = zeros(nvars);
nodesize = entropy{1}.H_x_2./log2(mi.N);

if time_pick ==1 && seg ==1 %plot dominant lags, all-segment averages

    I_normbyH = AllStats.I_normbyH;
    I_lag = AllStats.I_lag;
    TE = AllStats.TE;
    TI = AllStats.TI;
    U =  AllStats.U;
    S = AllStats.S;
    R = AllStats.R;
     
    if mi.nsegs>1
    H = AllStats.H_I;
    end
    
elseif seg==1 %plot specific lag, all-segment averages
      
    
    lag = time_pick-1; %if lag = 1 --> zero lag
    
    if lag>1
    I_lag = ones(nvars).*lagvect(lag-1);
    else %zero lag
    I_lag = zeros(nvars);
    end
    
    if mi.nsegs>1

        if lag>1
        I_normbyH = reshape(AllStats.I_normbyH_t(lag-1,:,:),nvars,nvars);
        else
          I_normbyH = AllStats.I_inst_normbyH;  
        end
 
        if lag>1
            TE = reshape(AllStats.TE_t(lag-1,:,:),nvars,nvars);
            TI = AllStats.TI;
            TI(AllStats.I_lag~=lagvect(lag-1))=0;
            U = AllStats.U;
            U(AllStats.I_lag~=lagvect(lag-1))=0;
            S = AllStats.S;
            S(AllStats.I_lag~=lagvect(lag-1))=0;
            R = AllStats.R;
            R(AllStats.I_lag~=lagvect(lag-1))=0;
        end
    else
        
        if lag>1
        I_normbyH = reshape(entropy{1}.I_normbyH(lag-1,:,:),nvars,nvars);
        else
        I_normbyH = entropy{1}.I_inst_normbyH;
        end
            
        if lag>1
            
            TE = reshape(entropy{1}.TE(lag-1,:,:),nvars,nvars);
            TI = entropy{1}.TI_T;
            TI(entropy{1}.I_dom_tau~=lagvect(lag-1))=0;
            U =  entropy{1}.U_T; 
            U(entropy{1}.I_dom_tau~=lagvect(lag-1))=0;
            S =  entropy{1}.S_T; 
            S(entropy{1}.I_dom_tau~=lagvect(lag-1))=0;
            R =  entropy{1}.R_T; 
            R(entropy{1}.I_dom_tau~=lagvect(lag-1))=0;
        end
    end
    
elseif time_pick == 1 && seg>1 %plot dominant lags, specific segment
    s = seg-1;
    I_normbyH = entropy{s}.I_dom;
    I_lag = entropy{s}.I_dom_tau;
    TE = entropy{s}.TE_T_normbyItot;
    TI = entropy{s}.TI_T;
    U =  entropy{s}.U_T_normbyItot;
    S =  entropy{s}.S_T_normbyItot;
    R =  entropy{s}.R_T_normbyItot;
    
else %plot specific lag, specific segment     
    lag = time_pick-1;
    s = seg-1;
  
    if lag>1
    I_normbyH = reshape(entropy{s}.I_normbyH(lag,:,:),nvars,nvars);
    I_lag = ones(nvars).*lagvect(lag-1);
    else
    I_normbyH = entropy{s}.I_inst_normbyH;
    I_lag = zeros(nvars); %zero-lag
    end
    
    if lag>1
    TE = entropy{s}.TE_T;
    TE(entropy{s}.I_dom_tau~=lagvect(lag-1))=0;
    TI = entropy{s}.TI_T;
    TI(entropy{s}.I_dom_tau~=lagvect(lag-1))=0;
    U = entropy{s}.U_T;
    U(entropy{s}.I_dom_tau~=lagvect(lag-1))=0;
    S = entropy{s}.S_T;
    S(entropy{s}.I_dom_tau~=lagvect(lag-1))=0;
    R = entropy{s}.R_T;
    R(entropy{s}.I_dom_tau~=lagvect(lag-1))=0;
    end
    
end

titletext = '';
uplim = .001;
bothways = 1;

if type==1       
        Links = I_normbyH;    
elseif type==2
        Links = TE; 
elseif type==3
        Links = TI;
elseif type==4
        Links = U;      
elseif type ==5
        Links = H; 
elseif type==6
        Links = S;
elseif type==7
        Links = R;
end




%only plot certain nodes acc to transnode and recnode
origlinks = Links;
if transnode >1 && recnode==1 %plot all links from certain transmitter
    Links =zeros(nvars);
    Links(transnode-1,:)=origlinks(transnode-1,:);
elseif transnode ==1 && recnode >1
    Links =zeros(nvars);
    Links(:,recnode-1)=origlinks(:,recnode-1);
elseif transnode>1 && recnode>1
    Links =zeros(nvars);
    Links(transnode-1,recnode-1)=origlinks(transnode-1,recnode-1);
end

for i =1:mi.nvars
    nodesize(i) = Links(i,i);
end

Weights = I_lag;
WeightRange = [0 max(lagvect) length(lagvect)+1];
A=1;
nodesize(nodesize<=0)=.01;
colorlabel = 'lag time';

nodesize(isnan(nodesize))=0.05;

fh=handles.plot_circle;
cla(fh);
PlotFunGUI(fh, Links,Weights, WeightRange, nodesize, varnames,titletext,A,uplim,bothways);

if mi.nsegs>1
axes(handles.plot_windows)

%show node pair relations if certain node pairs selected
contents=cellstr(get(handles.menu_transnode,'String'));
pairstring = contents{get(handles.menu_transnode,'Value')};

axesHandlesToChildObjects = findobj(gca, 'Type', 'line','Color','r');
	if ~isempty(axesHandlesToChildObjects)
		delete(axesHandlesToChildObjects);
    end
   
if seg>1   
    segline=line([seg-1 seg-1],[0 1],'LineWidth',3,'Color','r');
end

legend('Lagged I','Inst I','T_E','T/I','Location','East')

end


% --- Executes during object creation, after setting all properties.
function buttons_plottype_CreateFcn(hObject, eventdata, handles)
guidata(hObject,handles);


% --- Executes on button press in closepdf.
function closepdf_Callback(hObject, eventdata, handles)

close EntropyGUI_plot
wind = findobj(0,'Tag','EntropyGUI_main');
set(wind,'Visible','On')


% --- Executes on selection change in menu_transnode.
function menu_transnode_Callback(hObject, eventdata, handles)

set(handles.button_plot,'backgroundcolor','y')

% --- Executes during object creation, after setting all properties.
function menu_transnode_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menu_plottype.
function menu_plottype_Callback(hObject, eventdata, handles)

set(handles.button_plot,'backgroundcolor','y')

lag = get(handles.menu_lags,'Value');
type = get(handles.menu_plottype,'Value');
if lag>1 && type > 2 %T/I =3, U=4, H=5
    msgbox('Warning: T/I, U, H(Link) only for dominant time lags, reverting to dominant lags',...
        'error','error')
    set(handles.menu_lags,'Value',1)   
end

% --- Executes during object creation, after setting all properties.
function menu_plottype_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_save.
function button_save_Callback(hObject, eventdata, handles)

saveas(gca,'networkfig.jpg')


% --- Executes on selection change in menu_recnode.
function menu_recnode_Callback(hObject, eventdata, handles)

set(handles.button_plot,'backgroundcolor','y')

% --- Executes during object creation, after setting all properties.
function menu_recnode_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
