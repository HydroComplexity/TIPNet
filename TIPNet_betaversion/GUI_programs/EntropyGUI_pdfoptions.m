function varargout = EntropyGUI_pdfoptions(varargin)
% ENTROPYGUI_PDFOPTIONS MATLAB code for EntropyGUI_pdfoptions.fig
%      ENTROPYGUI_PDFOPTIONS, by itself, creates a new ENTROPYGUI_PDFOPTIONS or raises the existing
%      singleton*.
%
%      H = ENTROPYGUI_PDFOPTIONS returns the handle to a new ENTROPYGUI_PDFOPTIONS or the handle to
%      the existing singleton*.
%
%      ENTROPYGUI_PDFOPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENTROPYGUI_PDFOPTIONS.M with the given input arguments.
%
%      ENTROPYGUI_PDFOPTIONS('Property','Value',...) creates a new ENTROPYGUI_PDFOPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EntropyGUI_pdfoptions_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EntropyGUI_pdfoptions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EntropyGUI_pdfoptions

% Last Modified by GUIDE v2.5 09-May-2016 12:27:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EntropyGUI_pdfoptions_OpeningFcn, ...
                   'gui_OutputFcn',  @EntropyGUI_pdfoptions_OutputFcn, ...
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


% --- Executes just before EntropyGUI_pdfoptions is made visible.
function EntropyGUI_pdfoptions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EntropyGUI_pdfoptions (see VARARGIN)

% Choose default command line output for EntropyGUI_pdfoptions
handles.output = hObject;

handles.seg=1;

load('GUI_programs/Temps/projectname.mat')
load(fullpathname) %project file: contains mi structure of model information, data_orig, data

set(hObject,'toolbar','figure');

%initialize list boxes with node choices for plotting
%load('GUI_programs/Temps/dataset.mat')
choices1 = varnames; 
choices2 = ['-- (1D pdf)',varnames]; 
choices3 = ['-- (1D or 2D pdf)',varnames];
set(handles.choose_node1,'String',choices1)
set(handles.choose_node2,'String',choices2)
set(handles.choose_node3,'String',choices3)

%initialize list boxes for lags on second two node options
choices = ['0 (no lag)',num2cell(mi.lagvect)];
set(handles.choose_lag1,'String',choices);
set(handles.choose_lag2,'String',choices);

set(handles.choose_lag1,'Value',1)
set(handles.choose_lag2,'Value',1)


set(handles.choose_node1,'Value',1)
set(handles.choose_node2,'Value',1)
set(handles.choose_node3,'Value',1)


if strcmp(mi.bin_scheme,'global')==1
set(handles.local_global,'Value',1);
else
set(handles.local_global,'Value',2);
end

if strcmp(mi.method,'KDE')==1 %KDE method choses
set(handles.choose_method,'Value',1);
set(handles.hpanel,'Visible','on');
else %default to fixed bin method
set(handles.choose_method,'Value',2);
set(handles.hpanel,'Visible','off');   
end

%show h (smoothing parameter)for first node, hide other slider boxes

set(handles.show_h1,'String',num2str(mi.KDEparams.h1(1)));

set(handles.slider_h1,'Min',mi.KDEparams.h1_seg_orig(1)./5);
set(handles.slider_h1,'Max',mi.KDEparams.h1_seg_orig(1).*5);
set(handles.slider_h1,'Value',mi.KDEparams.h1(1));

if strcmp(mi.method,'fixed')==1
    set(handles.hpanel,'Visible','off')
    set(handles.slider_h1,'Visible','off');
%    set(handles.text_h1,'Visible','off');
end

set(handles.slider_h2,'Visible','Off');
set(handles.slider_h3,'Visible','Off');
set(handles.show_h2,'Visible','Off');
set(handles.text_h2,'Visible','Off');
set(handles.text_h3,'Visible','Off');

set(handles.setN,'String',num2str(mi.N))

nodes = [1 0 0];
save('GUI_programs/Temps/PDFnodes.mat','nodes')

%plot 1D pdf of first node, first segment
pdf = compute_pdfGUI(netdata{1}.data_process2(:,1),mi.N,mi.bin_scheme, mi.Range(:,1),mi.method,mi.KDEparams.h1(1));
axes(handles.plot_pdf)

info = compute_info_measures2(pdf);

bar(1/mi.N:1/mi.N:1,pdf);
xlim([-.01 1.1])

ylabel(sprintf('p(%s)',varnames{nodes(1)}))
xlabel(sprintf('%s',varnames{nodes(1)}))

%populate other boxes with entropy values
pdf = compute_pdfGUI(netdata{1}.data_process2(:,2),mi.N,mi.bin_scheme, mi.Range(:,2),mi.method,mi.KDEparams.h1(2));
info = compute_info_measures2(pdf);


%set(handles.textH2,'String',sprintf('H(%s)=%5.4f',varnames{2},info.Hx));

if mi.nvars>2
pdf = compute_pdfGUI(netdata{1}.data_process2(:,3),mi.N,mi.bin_scheme, mi.Range(:,3),mi.method,mi.KDEparams.h1(3));
info = compute_info_measures2(pdf);
%set(handles.textH3,'String',sprintf('H(%s)=%5.4f',varnames{3},info.Hx));
end
%set segment to first segment, or make invisible if only 1 segment
if mi.nsegs==1
    set(handles.pick_seg,'Visible','off')
    set(handles.text16,'Visible','off')
else
    set(handles.pick_seg,'Visible','on')
    set(handles.text16,'Visible','on')
    choices=[];
    for i =1:mi.nsegs
    choices{i} = num2str(i);
    end
    set(handles.pick_seg,'String',choices)
    
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EntropyGUI_pdfoptions wait for user response (see UIRESUME)
% uiwait(handles.EntropyGUI_pdf);


% --- Outputs from this function are returned to the command line.
function varargout = EntropyGUI_pdfoptions_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


% --- Executes on selection change in choose_node1.
function choose_node1_Callback(hObject, eventdata, handles)

load('GUI_programs/Temps/PDFnodes.mat')
nodes(1) = get(hObject,'Value');
save('GUI_programs/Temps/PDFnodes.mat','nodes')

load('GUI_programs/Temps/projectname.mat')
load(fullpathname) %project file: contains mi structure of model information, data_orig, data


%set h(node 1) value to 1D, 2D, or 3D value
if nodes(2)==0
    h=mi.KDEparams.h1(nodes(1));
    set(handles.slider_h1,'Min',mi.KDEparams.h1_seg_orig(nodes(1))./5);
    set(handles.slider_h1,'Max',mi.KDEparams.h1_seg_orig(nodes(1)).*5);
elseif nodes(3)==0
    h=mi.KDEparams.h2(nodes(1));
    set(handles.slider_h1,'Min',mi.KDEparams.h2_seg_orig(nodes(1))./5);
    set(handles.slider_h1,'Max',mi.KDEparams.h2_seg_orig(nodes(1)).*5);
else
    h=mi.KDEparams.h3(nodes(1));
    set(handles.slider_h1,'Min',mi.KDEparams.h3_seg_orig(nodes(1))./5);
    set(handles.slider_h1,'Max',mi.KDEparams.h3_seg_orig(nodes(1)).*5);
end

set(handles.show_h1,'String',num2str(h));
set(handles.slider_h1,'Value',h);

% --- Executes during object creation, after setting all properties.
function choose_node1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in choose_node2.
function choose_node2_Callback(hObject, eventdata, handles)

load('GUI_programs/Temps/PDFnodes.mat')

load('GUI_programs/Temps/projectname.mat')
load(fullpathname) %project file: contains mi structure of model information, data_orig, data


nodes(2) = get(hObject,'Value')-1;
if nodes(2)==0 %1D pdf
    set(handles.choose_node3,'Value',1)
    nodes(3)=0;     
    set(handles.show_h2,'Visible','Off');
    set(handles.text_h2,'Visible','Off');
    set(handles.slider_h2,'Visible','Off');
    set(handles.show_h3,'Visible','Off');
    set(handles.text_h3,'Visible','Off');
    set(handles.slider_h3,'Visible','Off');
end

if nodes(2)==0    
    h1 = mi.KDEparams.h1(nodes(1));
    set(handles.slider_h1,'Min',mi.KDEparams.h1_seg_orig(nodes(1))./5);
    set(handles.slider_h1,'Max',mi.KDEparams.h1_seg_orig(nodes(1)).*5);
elseif nodes(3)==0  
    h1 = mi.KDEparams.h2(nodes(1));
    h2 = mi.KDEparams.h2(nodes(2));
    set(handles.slider_h1,'Min',mi.KDEparams.h2_seg_orig(nodes(1))./5);
    set(handles.slider_h1,'Max',mi.KDEparams.h2_seg_orig(nodes(1)).*5);
    set(handles.slider_h2,'Min',mi.KDEparams.h2_seg_orig(nodes(2))./5);
    set(handles.slider_h2,'Max',mi.KDEparams.h2_seg_orig(nodes(2)).*5);
else
    h1 = mi.KDEparams.h3(nodes(1));
    h2 = mi.KDEparams.h3(nodes(2));
    h3 = mi.KDEparams.h3(nodes(3));
    set(handles.slider_h1,'Min',mi.KDEparams.h3_seg_orig(nodes(1))./5);
    set(handles.slider_h1,'Max',mi.KDEparams.h3_seg_orig(nodes(1)).*5);
    set(handles.slider_h2,'Min',mi.KDEparams.h3_seg_orig(nodes(2))./5);
    set(handles.slider_h2,'Max',mi.KDEparams.h3_seg_orig(nodes(2)).*5);
    set(handles.slider_h3,'Min',mi.KDEparams.h3_seg_orig(nodes(3))./5);
    set(handles.slider_h3,'Max',mi.KDEparams.h3_seg_orig(nodes(3)).*5);    
end
    
    %set h(node 1) value to 1D value
    set(handles.show_h1,'String',num2str(h1));
    set(handles.slider_h1,'Value',h1);
    
if nodes(2)>0 %2D or 3D pdf
    set(handles.text_h2,'Visible','On');
    set(handles.show_h2,'Visible','On');
    set(handles.slider_h2,'Visible','On');
   
    set(handles.show_h2,'String',num2str(h2));
    set(handles.slider_h2,'Value',h2);  
end 

if nodes(3)>0 %2D or 3D pdf
    set(handles.text_h3,'Visible','On');
    set(handles.show_h3,'Visible','On');
    set(handles.slider_h3,'Visible','On');
   
    set(handles.show_h3,'String',num2str(h3));
    set(handles.slider_h3,'Value',h3);  
end 

save('GUI_programs/Temps/PDFnodes.mat','nodes')

% --- Executes during object creation, after setting all properties.
function choose_node2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in choose_node3.
function choose_node3_Callback(hObject, eventdata, handles)

load('GUI_programs/Temps/PDFnodes.mat')

load('GUI_programs/Temps/projectname.mat')
load(fullpathname) %project file: contains mi structure of model information, data_orig, data
contents = cellstr(get(hObject,'String'));
nodes(3) = get(hObject,'Value')-1;

save('GUI_programs/Temps/PDFnodes.mat','nodes')

if nodes(3)==0 %1D or 2D pdf
    set(handles.show_h3,'Visible','Off');
    set(handles.text_h3,'Visible','Off');
    set(handles.slider_h3,'Visible','Off');
end

if nodes(3)==0  && nodes(2)==0  
    h1 = mi.KDEparams.h1(nodes(1));
    set(handles.slider_h1,'Min',mi.KDEparams.h1_seg_orig(nodes(1))./5);
    set(handles.slider_h1,'Max',mi.KDEparams.h1_seg_orig(nodes(1)).*5);
 
elseif nodes(3)==0  
    h1 = mi.KDEparams.h2(nodes(1));
    h2 = mi.KDEparams.h2(nodes(2));
    
    set(handles.slider_h1,'Min',mi.KDEparams.h2_seg_orig(nodes(1))./5);
    set(handles.slider_h1,'Max',mi.KDEparams.h2_seg_orig(nodes(1)).*5);
    set(handles.slider_h2,'Min',mi.KDEparams.h2_seg_orig(nodes(2))./5);
    set(handles.slider_h2,'Max',mi.KDEparams.h2_seg_orig(nodes(2)).*5);
   
else
    h1 = mi.KDEparams.h3(nodes(1));
    h2 = mi.KDEparams.h3(nodes(2));
    h3 = mi.KDEparams.h3(nodes(3));
    
    set(handles.slider_h1,'Min',mi.KDEparams.h3_seg_orig(nodes(1))./5);
    set(handles.slider_h1,'Max',mi.KDEparams.h3_seg_orig(nodes(1)).*5);
    set(handles.slider_h2,'Min',mi.KDEparams.h3_seg_orig(nodes(2))./5);
    set(handles.slider_h2,'Max',mi.KDEparams.h3_seg_orig(nodes(2)).*5);
    set(handles.slider_h3,'Min',mi.KDEparams.h3_seg_orig(nodes(3))./5);
    set(handles.slider_h3,'Max',mi.KDEparams.h3_seg_orig(nodes(3)).*5);   
end
    
    %set h(node 1) value to 1D value
    set(handles.show_h1,'String',num2str(h1));
    set(handles.slider_h1,'Value',h1);
    
if nodes(2)>0 %2D or 3D pdf
    set(handles.text_h2,'Visible','On');
    set(handles.show_h2,'Visible','On');
    set(handles.slider_h2,'Visible','On');
   
    set(handles.show_h2,'String',num2str(h2));
    set(handles.slider_h2,'Value',h2);  
end 

if nodes(3)>0 %3D pdf
    set(handles.text_h3,'Visible','On');
    set(handles.show_h3,'Visible','On');
    set(handles.slider_h3,'Visible','On');
    set(handles.show_h3,'String',num2str(h3));
    set(handles.slider_h3,'Value',h3);  
end 

% --- Executes during object creation, after setting all properties.
function choose_node3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pdfbutton.
function pdfbutton_Callback(hObject, eventdata, handles)

clear handles.Colorbar


load('GUI_programs/Temps/PDFnodes.mat')
load('GUI_programs/Temps/projectname.mat')
load(fullpathname) %project file: contains mi structure of model information, data_orig, data

cla(handles.plot_pdf,'reset')

seg=handles.seg;
method = mi.method;
N=mi.N;

%lag times selected for nodes 2 and 3
l1 = get(handles.choose_lag1,'Value');
if l1 ==1
    lag1 =0;
else
    lag1 = mi.lagvect(l1-1);
end
l2 = get(handles.choose_lag2,'Value');
if l2 ==1
    lag2 =0;
else
    lag2 = mi.lagvect(l2-1);
end


if nodes(2)==0 %compute and plot 1D pdf (line plot)
       
    dat = netdata{seg}.data_process2(:,nodes(1));
    r = mi.Range(:,nodes(1));
    h = mi.KDEparams.h1(nodes(1));
    pdf = compute_pdfGUI(dat,N,mi.bin_scheme, r,method,h);
    axes(handles.plot_pdf)
    bar(1/N:1/N:1,pdf)
    xlim([-.01 1.1])

    ylabel(sprintf('p(%s)',varnames{nodes(1)}))
    xlabel(sprintf('%s',varnames{nodes(1)}))
    
    info = compute_info_measures2(pdf);

    
elseif nodes(3)==0 %compute and plot 2D pdf (surface plot)
   
    dat1 = netdata{seg}.data_process2((lag1+1):end,nodes(1));
    len = length(dat1);
    dat2 = netdata{seg}.data_process2(1:len,nodes(2)); %lagged dat2 by lag1
    
    dat = [dat1 dat2];
    r = [mi.Range(:,nodes(1)) mi.Range(:,nodes(2))];
    h = [mi.KDEparams.h2(nodes(1)) mi.KDEparams.h2(nodes(2))];
    pdf = compute_pdfGUI(dat,mi.N,mi.bin_scheme, r,method,h);
    axes(handles.plot_pdf)
    
    pdf(pdf==0)=nan;
    im = imagesc(1/N:1/N:1,1/N:1/N:1,pdf);
    set(im, 'AlphaData', ~isnan(pdf))
    
    set(gca,'Ydir','normal');
    axis square
    colorbar

    xlabel(sprintf('%s(t-%d)',varnames{nodes(2)},lag1))
    ylabel(sprintf('%s',varnames{nodes(1)}))
    title(sprintf('p(%s, %s)',varnames{nodes(1)},varnames{nodes(2)}))
    
else %compute and plot 3D pdf (cloud plot)
    
    
    if lag2>lag1
         
         dat1 = netdata{seg}.data_process2((lag2+1):end,nodes(1));
         len = length(dat1);
         dat2 = netdata{seg}.data_process2((lag2-lag1+1):(lag2-lag1+len),nodes(2)); %lagged dat2 by lag1
         dat3 = netdata{seg}.data_process2(1:len,nodes(3));
    elseif lag1>lag2
         dat1 = netdata{seg}.data_process2((lag1+1):end,nodes(1));
         len = length(dat1);
         dat2 = netdata{seg}.data_process2((lag1-lag2+1):(lag1-lag2+len),nodes(2)); %lagged dat2 by lag1
         dat3 = netdata{seg}.data_process2(1:len,nodes(3));
    else
         dat1 = netdata{seg}.data_process2((lag1+1):end,nodes(1));
         len = length(dat1);
         dat2 = netdata{seg}.data_process2(1:len,nodes(2)); %lagged dat2 by lag1
         dat3 = netdata{seg}.data_process2(1:len,nodes(3));
    end
    
    dat = [dat1 dat2 dat3];
    r = [mi.Range(:,nodes(1)) mi.Range(:,nodes(2)) mi.Range(:,nodes(3))];
    h = [mi.KDEparams.h3(nodes(1)) mi.KDEparams.h3(nodes(2)) mi.KDEparams.h3(nodes(3))];
    pdf = compute_pdfGUI(dat,mi.N,mi.bin_scheme, r,method,h); 
    
    N = size(pdf,1);
    
    step=1;
    if N>15 && N<30
    step = 2;
    elseif N>30
    step = 3;
    end
    cvect = jet(100);
    
    minpdf = min(min(min(pdf(pdf~=0))));
    maxpdf = max(max(max(pdf)));

    axes(handles.plot_pdf)
    cla(gca,'reset')
    hold on;
    grid on
    
     labels = num2cell(round(1000.*linspace(minpdf,maxpdf,8))./1000);
     colormap(jet(8));
     caxis([minpdf maxpdf]);
     c=colorbar;
    
    for ii = 1:step:N
        ci = (ii-1)/(N-1);
        for jj=1:step:N
            cj = (jj-1)/(N-1);
            for kk = 1:step:N
                ck = (kk-1)/(N-1);
                
                pdfcat=(pdf(ii,jj,kk)-minpdf)./(maxpdf-minpdf);
                pdfcat = round(pdfcat.*100);
                if pdfcat > 1
                    plot3(ci,cj,ck,'o','Color',cvect(pdfcat,:),'LineWidth',pdfcat./30);
                end
                
            end
        end
    end
 
    xlabel(sprintf('%s',varnames{nodes(1)})); 
    ylabel(sprintf('%s(t-%d)',varnames{nodes(2)},lag1));
    zlabel(sprintf('%s(t-%d)',varnames{nodes(3)},lag2));
  
    hold off;

axis square
xlim([0 1]); ylim([0 1]); zlim([0 1]);
view(20,10)   
end

%for 2D or 3D pdf, compute mutual information and statisitical significance
% if nodes(2)>0
%     dat = [netdata{1}.data_process2(:,nodes(1)) netdata{1}.data_process2(:,nodes(2))];
%     r = [mi.Range(:,nodes(1)) mi.Range(:,nodes(2))];
%     h = [mi.KDEparams.h2(nodes(1)) mi.KDEparams.h2(nodes(2))];
%     pdf = compute_pdfGUI(dat,mi.N,mi.bin_scheme, r,method,h);
%     info = compute_info_measures2(pdf);
%     
%     I_lagged = info.I;
%     
%     I_shuff=zeros(1,10);
%     for test=1:10
%         shuff1 = randsample(netdata{1}.data_process2(:,nodes(1)),length(netdata{1}.data_process2(:,nodes(1))));
%         shuff2 = randsample(netdata{1}.data_process2(:,nodes(2)),length(netdata{1}.data_process2(:,nodes(2))));
%         dat = [shuff1 shuff2];
%         r = [mi.Range(:,nodes(1)) mi.Range(:,nodes(2))];
%         h = [mi.KDEparams.h2(nodes(1)) mi.KDEparams.h2(nodes(2))];
%         pdf = compute_pdfGUI(dat,mi.N,mi.bin_scheme, r,method,h);
%         infoshuff = compute_info_measures2(pdf);
%         I_shuff(test) = infoshuff.I;
%     end
%     I_shuff_sig = mean(I_shuff)+3*std(I_shuff); 
%     if I_shuff_sig>I_lagged
%         I_lagged=0;
%     end
%     %I_lagged = max(0,I_lagged-I_shuff_sig);
%     
%     %if significant I_lagged, and 3d pdf chosen, compute other measures
%     if  nodes(3)>0
%         %Xt (node1), Yw (node3), Yf (node2)
%         dat = [netdata{1}.data_process2(:,nodes(1)) netdata{1}.data_process2(:,nodes(3)) netdata{1}.data_process2(:,nodes(2))];
%         r = [mi.Range(:,nodes(1)) mi.Range(:,nodes(3)) mi.Range(:,nodes(2))];
%         h = [mi.KDEparams.h3(nodes(1))  mi.KDEparams.h3(nodes(3)) mi.KDEparams.h3(nodes(2))];
%         pdf = compute_pdfGUI(dat,mi.N,mi.bin_scheme, r,method,h);
%         info3d = compute_info_measures2(pdf);
%         
%    
% 
%     end

  
  %save pdf in Temps
  save('GUI_programs/Temps/pdf.mat','pdf')
        
%end


function setN_Callback(hObject, eventdata, handles)

N=str2double(get(hObject,'String'));

if N>100
h = msgbox('max N=100', 'Error','error'); 
N=100;
elseif N<2
h = msgbox('min N=2', 'Error','error'); 
N=2;
elseif mod(N,1)>0
N=round(N);
end

load('GUI_programs/Temps/projectname.mat')
load(fullpathname)
mi.N = N;
save(fullpathname,'mi','-append')

set(hObject,'String',num2str(N))

% --- Executes during object creation, after setting all properties.
function setN_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in local_global.
function local_global_Callback(hObject, eventdata, handles)

load('GUI_programs/Temps/projectname.mat')
load(fullpathname)

index_selected = get(hObject,'Value');

if index_selected==1
    mi.bin_scheme = 'global';
else
    mi.bin_scheme = 'local';
end

save(fullpathname,'mi','-append')

% --- Executes during object creation, after setting all properties.
function local_global_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_h1_Callback(hObject, eventdata, handles)

h = get(hObject,'Value');
set(handles.show_h1,'String',num2str(h))

load('GUI_programs/Temps/projectname.mat')
load(fullpathname)

load('GUI_programs/Temps/PDFnodes.mat')

if nodes(3)==0 && nodes(2)==0 %1D pdf
    mi.KDEparams.h1(nodes(1))=h;
elseif nodes(3)==0 %2D pdf
    mi.KDEparams.h2(nodes(1))=h;
else %3D pdf
    mi.KDEparams.h3(nodes(1))=h;
end

save(fullpathname,'mi','-append')

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider_h1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_h1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_h2_Callback(hObject, eventdata, handles)
% hObject    handle to slider_h2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
h = get(hObject,'Value');
set(handles.show_h2,'String',num2str(h))

load('GUI_programs/Temps/projectname.mat')
load(fullpathname)

load('GUI_programs/Temps/PDFnodes.mat')

if nodes(3)==0 %2D pdf
    mi.KDEparams.h2(nodes(2))=h;
else %3D pdf
    mi.KDEparams.h3(nodes(2)) =h;
end

save(fullpathname,'mi','-append')

% --- Executes during object creation, after setting all properties.
function slider_h2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_h2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_h3_Callback(hObject, eventdata, handles)
% hObject    handle to slider_h3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

h = get(hObject,'Value');
set(handles.show_h3,'String',num2str(h))

load('GUI_programs/Temps/projectname.mat')
load(fullpathname)

load('GUI_programs/Temps/PDFnodes.mat')

mi.KDEparams.h3(nodes(3))=h;

save(fullpathname,'mi','-append')

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider_h3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_h3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in button_reset.
function button_reset_Callback(hObject, eventdata, handles)
% reset h (smoothing params) back to original values
load('GUI_programs/Temps/projectname.mat')
load('GUI_programs/Temps/PDFnodes.mat')
load(fullpathname)

mi.KDEparams.h1=mi.KDEparams.h1_seg_orig;
mi.KDEparams.h2=mi.KDEparams.h2_seg_orig;
mi.KDEparams.h3=mi.KDEparams.h3_seg_orig;

if nodes(2)==0 && nodes(3)==0 %1D pdf
h1 = mi.KDEparams.h1(nodes(1));
set(handles.show_h1,'String',num2str(h1));
    set(handles.slider_h1,'Value',h1);
    
elseif nodes(3)==0 %2D pdf
    h1 = mi.KDEparams.h2(nodes(1));
    h2 = mi.KDEparams.h2(nodes(2));
    set(handles.show_h1,'String',num2str(h1));
    set(handles.slider_h1,'Value',h1);
    set(handles.show_h2,'String',num2str(h2));
    set(handles.slider_h2,'Value',h2);
else %3D pdf
    h1 = mi.KDEparams.h3(nodes(1));
    h2 = mi.KDEparams.h3(nodes(2));
    h3 = mi.KDEparams.h3(nodes(3));
    set(handles.show_h1,'String',num2str(h1));
    set(handles.slider_h1,'Value',h1);
    set(handles.show_h2,'String',num2str(h2));
    set(handles.slider_h2,'Value',h2);
    set(handles.show_h3,'String',num2str(h3));
    set(handles.slider_h3,'Value',h3);
end


save(fullpathname,'mi','-append');


% --- Executes on button press in button_done.
function button_done_Callback(hObject, eventdata, handles)
% hObject    handle to button_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close EntropyGUI_pdfoptions

wind = findobj(0,'Tag','EntropyGUI_main');
set(wind,'Visible','On')


% --- Executes on selection change in pick_seg.
function pick_seg_Callback(hObject, eventdata, handles)

seg=get(hObject,'Value');

handles.seg=seg;

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function pick_seg_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in choose_method.
function choose_method_Callback(hObject, eventdata, handles)

load('GUI_programs/Temps/projectname.mat')
load(fullpathname)

index_selected = get(hObject,'Value');

if index_selected==1
    mi.method = 'KDE';
    set(handles.hpanel,'Visible','on')
    set(handles.slider_h1,'Visible','on');
    
    
else
    mi.method = 'fixed';
    set(handles.hpanel,'Visible','off')
end

save(fullpathname,'mi','-append')

% --- Executes during object creation, after setting all properties.
function choose_method_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_savefig.
function button_savefig_Callback(hObject, eventdata, handles)
% hObject    handle to button_savefig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Save pdf to UserData folder
load('GUI_programs/Temps/pdf.mat');
load('GUI_programs/Temps/pdfnodes.mat');
save('UserData/pdf_saved.mat','pdf','nodes');

% --- Executes on selection change in choose_lag1.
function choose_lag1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function choose_lag1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in choose_lag2.
function choose_lag2_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function choose_lag2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
