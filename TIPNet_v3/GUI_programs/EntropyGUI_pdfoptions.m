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
else %default to fixed bin method
set(handles.choose_method,'Value',2);  
end

set(handles.setN,'String',num2str(mi.N))

nodes = [1 0 0];
save('GUI_programs/Temps/PDFnodes.mat','nodes')

z_eff = mi.DataPrep.Z_effect;

%plot 1D pdf of first node, first segment
pdf = compute_pdfGUI(netdata{1}.data_process2(:,1),mi.N,mi.bin_scheme, mi.Range(:,1),mi.method,z_eff(1+1));
axes(handles.plot_pdf)

info = compute_info_measures(pdf);

bar(1/mi.N:1/mi.N:1,pdf);
xlim([-.01 1.1])

ylabel(sprintf('p(%s)',varnames{nodes(1)}))
xlabel(sprintf('%s',varnames{nodes(1)}))

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
z_eff = mi.DataPrep.Z_effect;

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
    z = z_eff(nodes(1)+1);
    pdf = compute_pdfGUI(dat,N,mi.bin_scheme, r,method,z);
    axes(handles.plot_pdf)
    bar(1/N:1/N:1,pdf)
    xlim([-.01 1.1])

    ylabel(sprintf('p(%s)',varnames{nodes(1)}))
    xlabel(sprintf('%s',varnames{nodes(1)}))
    
    info = compute_info_measures(pdf);
    
elseif nodes(3)==0 %compute and plot 2D pdf (surface plot)
    
    
    dat1 = netdata{seg}.data_process2((lag1+1):end,nodes(1));
    len = length(dat1);
    dat2 = netdata{seg}.data_process2(1:len,nodes(2)); %lagged dat2 by lag1
    
    dat = [dat1 dat2];
    r = [mi.Range(:,nodes(1)) mi.Range(:,nodes(2))];
    z = [z_eff(nodes(1)+1) z_eff(nodes(2)+1)];
    pdf = compute_pdfGUI(dat,mi.N,mi.bin_scheme, r,method,z);
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
    z = [z_eff(nodes(1)+1) z_eff(nodes(2)+1) z_eff(nodes(3)+1)];
    pdf = compute_pdfGUI(dat,mi.N,mi.bin_scheme, r,method,z); 
    
    N = size(pdf,1);
    
    step=1;

    cvect = jet(1000);
    
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
                pdfcat = round(pdfcat.*1000);
                if pdfcat > 1
                    plot3(ci,cj,ck,'o','Color',cvect(pdfcat,:),'LineWidth',pdfcat./300);
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
else
    mi.method = 'fixed';
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

saveas(handles.plot_pdf,'pdffig.jpg')


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
