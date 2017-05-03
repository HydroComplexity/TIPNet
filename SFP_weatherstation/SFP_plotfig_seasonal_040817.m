
%SFP network plotting: seasonal trend for summer 2015

%top: plot of U, R, S total information over season (smoothed)
%bottom: plot of U,R,S from source pairs {Ta, Other} to RH, where Other =
%LWet, RH, WS, Rg, etc etc. (also smoothed)

clear all
close all
clc

filename = 'SFPproject_032717.mat';
load(filename);

PlotTemp =1;

nvars = length(varnames);

Yr = 2015;
start_DOY = 120;
end_DOY = 299;

for s =1:length(entropy)
   day(s)=floor(mean(netdata{s}.decdoy));
   year(s)=mean(floor(netdata{s}.decyear));
end

start_seg = find(day==start_DOY & year==Yr);
end_seg = find(day==end_DOY & year==Yr);
start_seg=start_seg(1);
end_seg = end_seg(end);

maxday = max(netdata{end_seg}.decdoy);
minday = min(netdata{start_seg}.decdoy);

Utotal = zeros(end_seg-start_seg+1,6);
Rtotal = zeros(end_seg-start_seg+1,6);
Stotal = zeros(end_seg-start_seg+1,6);

%% plot seasonal shifts in links

pairs = combnk(1:5,2);

ct=1;
cvect = jet(150);

cvect1 = autumn(120);
cvect2 = winter(120);

var =1; %target = RH
i =3; %source = Ta


dayct=1;

for s = start_seg:end_seg
             
    decyeartot(s) = mean(netdata{s}.decyear);
    PPTtot(s) = sum(netdata{s}.data_orig(:,5));
    decdoytot(s) = floor(mean(netdata{s}.decdoy));
    decdoy(s) = mean(netdata{s}.decdoy);
    
    if s==1
        DOY(dayct)= decdoytot(s);
        PPT(dayct) = PPTtot(s);
    elseif s>1 && decdoytot(s)==decdoytot(s-1)
    PPT(dayct)=PPT(dayct)+PPTtot(s);
    
    else
        
    dayct = dayct+1;
    DOY(dayct) = decdoytot(s);
    PPT(dayct)=PPTtot(s);
    end       
   
end


for s = start_seg:end_seg
    


    decyear(ct) = mean(netdata{s}.decyear);
    decdoy(ct) = max(netdata{s}.decdoy);

    StargetVar = entropy{s}.U(i,var);

    S1lag = entropy{s}.I_dom_lag(i,var);
    
    Utotal(s)=StargetVar;

end



%% synergistic and redundant info
    
for s = start_seg:end_seg
    decyear(ct) = mean(netdata{s}.decyear);
    decdoy(ct) = max(netdata{s}.decdoy);

    StargetVar = entropy{s}.S(i,var);
    StargetPair = entropy{s}.S_pair(i,var);
    RtargetVar = entropy{s}.R(i,var);
    RtargetPair = entropy{s}.R_pair(i,var);
  
    if StargetVar>0
    Stotal(s,StargetPair)=StargetVar;
    end   
    
    if RtargetVar>0
    Rtotal(s,RtargetPair)=RtargetVar;
    end
    ct=ct+1;
      
end
    
for s = start_seg:end_seg
    Stotal_all(s) = sum(sum(entropy{s}.S));
    Rtotal_all(s) = sum(sum(entropy{s}.R));
    Utotal_all(s) = sum(sum(entropy{s}.U));
end


%%

s = 4*15;

Utotalsum = cumsum(Utotal);
Rtotalsum = cumsum(Rtotal);
Stotalsum = cumsum(Stotal);

S_RH_sm = smooth(Stotal(:,1),s);
S_WS_sm = smooth(Stotal(:,2),s);
S_Rg_sm = smooth(Stotal(:,4),s);
S_LW_sm = smooth(Stotal(:,6),s);

R_RH_sm = smooth(Rtotal(:,1),s);
R_WS_sm = smooth(Rtotal(:,2),s);
R_Rg_sm = smooth(Rtotal(:,4),s);
R_LW_sm = smooth(Rtotal(:,6),s);

U_sm = smooth(Utotal(:,1),s);

z = zeros(size(S_RH_sm));
DOY = [decdoy, fliplr(decdoy)];

fig = figure(12);
h1 = axes('Units','centimeters','Position',[1 1 7 5]);
hold on

A = [z; flipud(S_RH_sm)];
B = [z; flipud(S_RH_sm+S_WS_sm)];
C = [z; flipud(S_RH_sm+S_WS_sm+ S_Rg_sm)];
D = [z; flipud(S_RH_sm+S_WS_sm+ S_Rg_sm+S_LW_sm)]; 
E = [z; flipud(S_RH_sm+S_WS_sm+ S_Rg_sm+S_LW_sm+R_RH_sm)]; 
F = [z; flipud(S_RH_sm+S_WS_sm+ S_Rg_sm+S_LW_sm+R_RH_sm+R_WS_sm)]; 
G = [z; flipud(S_RH_sm+S_WS_sm+ S_Rg_sm+S_LW_sm+R_RH_sm+R_WS_sm+R_Rg_sm)];
H = [z; flipud(S_RH_sm+S_WS_sm+ S_Rg_sm+S_LW_sm+R_RH_sm+R_WS_sm+R_Rg_sm+R_LW_sm)];
I = [z; flipud(S_RH_sm+S_WS_sm+ S_Rg_sm+S_LW_sm+R_RH_sm+R_WS_sm+R_Rg_sm+R_LW_sm+U_sm)];

cvect = [130 150 250; 100 190 130; 240 170 130; 230 130 140;...
60 80 160; 20 150 70; 250 130 30; 240 30 30; 30 40 90]./255;


fill(DOY, I,cvect(9,:))
fill(DOY, H,cvect(8,:))
fill(DOY, G,cvect(7,:))
fill(DOY, F,cvect(6,:))
fill(DOY, E,cvect(5,:))
fill(DOY, D,cvect(4,:))
fill(DOY, C,cvect(3,:))
fill(DOY, B,cvect(2,:))
fill(DOY, A,cvect(1,:))
xlim([minday+5 maxday-5])

h2 = axes('Units','centimeters','Position',[1 8 7 2.5]);
hold on


TotalI = Stotal_all+Rtotal_all+Utotal_all;

S_sm = smooth(Stotal_all,s);
R_sm = smooth(Rtotal_all,s);
U_sm = smooth(Utotal_all,s);


R = [z; flipud(R_sm+U_sm+S_sm)];
U = [z; flipud(S_sm+U_sm)];
S = [z; flipud(S_sm)];


fill(DOY,R,[0 0 0])
fill(DOY,U,[.5 .5 .5])
fill(DOY,S,[1 1 1])

xlim([minday+5 maxday-5])

