
%SFP network plotting from May/June 2016 SFP runs
%4 networks per day based on solar radiation (2 day and 2 night)
%goes from June 2014 to May 2016, 7 variables at 1 min resolution
%netdata{seg}.day_index variable = 1 for day time networks

clear all
close all
clc

filename = 'SFPproject_032717.mat';
load(filename);

for i=1:length(netdata)
    decdoy(i) = floor(mean(netdata{i}.decdoy));
end

PlotDOY = 124; %in Goodwell and Kumar 2017: 
%Plotting DOY 124 for LWET, 189, 197, and 293 for RH
ind = find(decdoy==PlotDOY);
start_seg = min(ind);
end_seg = max(ind);

maxday = max(netdata{end_seg}.decdoy);
minday = min(netdata{start_seg}.decdoy);


%% Point 1: plot 3-day timeseries, orginal data + filtered/normalized for network segments

for s =start_seg:end_seg

day = netdata{s}.decdoy;
origdata = netdata{s}.data_orig;
networkdata = netdata{s}.data_processed;

varvect = 1:6;
cvector = [60 80 160; 0 104 56; 100 45 145; 190 30 30; 240 90 40; 140 200 60]./255;
ct=1;
for i =1:6

figure(1)
subplot(6,2,ct)
hold on
plot(day,origdata(:,varvect(i)),'Color',cvector(i,:))
xlim([minday maxday])

subplot(6,2,ct+1)
hold on
plot(day,networkdata(:,varvect(i)),'Color',cvector(i,:))
%line([max(day) max(day)],[0 1],'Color','k')
xlim([minday maxday])
ct=ct+2;

end

end

figure(1)
ct=1;
for i =1:6
subplot(6,2,ct)
title(sprintf('Weather Station Data: %s',varnames{varvect(i)}))
subplot(6,2,ct+1)
title(sprintf('Pre-processed Data: %s',varnames{varvect(i)}))
ct=ct+2;
end

subplot(6,2,1)
ylim([.85 1.05])
ylabel('fraction')
subplot(6,2,3)
ylim([0 11])
ylabel('m/s')
subplot(6,2,5)
ylim([5 25])
ylabel('^{\circ} C')
subplot(6,2,7)
ylim([0 1200])
ylabel('W/m^2')
subplot(6,2,9)
ylim([0 1.2])
ylabel('mm/min')
subplot(6,2,11)
ylim([400 1100])
ylabel('counts')

%% network temporal plots: unique information
ct=1;
cvect = flipud(jet(61));
plotvar = [1 6];
%cvectlag = gray(10);

for v = 1:length(plotvar)
ct=1;
var=(plotvar(v));
figure(var+10)
subplot(3,1,1)
hold on

for s = start_seg:end_seg
    
    if var==1
    decyeartot(s) = mean(netdata{s}.decyear);
    PPTtot(s) = sum(netdata{s}.data_orig(:,6));
    end
    
    H = entropy{s}.H_x_2(var);
 
    if H >.005

    decyear(ct) = mean(netdata{s}.decyear);
    decdoy(ct) = max(netdata{s}.decdoy);

    for i =1:6
    offset = i./50;
    StargetVar = entropy{s}.U(i,var);
    if StargetVar>1
        Wid = 12;
    elseif StargetVar>.75
        Wid = 8;
    elseif StargetVar>.25
        Wid = 5;
    else
        Wid = 2;
    end
    
    %Wid = sqrt(StargetVar.*100);
    
    StargetVar_frac = entropy{s}.U_normbyItot(i,var);
    Wid2 = StargetVar_frac.*10;
    color2 = round(StargetVar_frac*100)+1;
    
   
    if StargetVar>0 %&& StargetVar==max(entropy{s}.U(:,var))
    S1lag = entropy{s}.I_dom_lag(i,var);
    
    figure(var+10)
    subplot(3,1,1)
    plot(decdoy(ct)-3/50, i,'.','Color',cvect(S1lag,:),'MarkerSize',Wid*5)
    
    end
    
    end
    
    
    
    ct=ct+1;
    end
    
    line([max(netdata{s}.decdoy) max(netdata{s}.decdoy)],[0 8],'Color','k')
 
end

figure(var+10)
subplot(3,1,1)
xlim([minday maxday])
ylim([.5 7.2])
set(gca,'Ytick',1:6,'Yticklabel',varnames)
title(sprintf('Unique Sources to %s',varnames{var}))
ylabel('dominant node pairs')


end

%% synergistic info
for v = 1:length(plotvar)
ct=1;

var=(plotvar(v));

figure(var+10)
subplot(3,1,2)
hold on

for s = start_seg:end_seg
      
    H = entropy{s}.H_x_2(var);

 
    if H >.005
    decyear(ct) = mean(netdata{s}.decyear);
    decdoy(ct) = max(netdata{s}.decdoy);

    for i =1:6
    offset = i./50;
    StargetVar = entropy{s}.S(i,var);
    if StargetVar>1
        Wid = 12;
    elseif StargetVar>.75
        Wid = 8;
    elseif StargetVar>.25
        Wid = 5;
    else
        Wid = 1;
    end
    
    StargetVar_frac = entropy{s}.S_normbyItot(i,var);
    Wid2 = StargetVar_frac.*10;
    color2 = round(StargetVar_frac*100)+1;
    
    
    if StargetVar>0 %&& StargetVar==max(entropy{s}.U(:,var))
    S1lag = entropy{s}.I_dom_lag(i,var);
    StargetPair = entropy{s}.S_pair(i,var);
    S2lag = entropy{s}.I_dom_lag(StargetPair,var);
    
    figure(var+10)
    subplot(3,1,2)
    line([decdoy(ct)-offset decdoy(ct)-offset],[i StargetPair],'LineWidth',Wid,'Color','k')
    plot(decdoy(ct)-offset, i,'.','Color',cvect(S1lag,:),'MarkerSize',Wid*5)
    plot(decdoy(ct)-offset,StargetPair,'.','Color',cvect(S2lag,:),'MarkerSize',Wid*5)
    end
    
    end
    
    
    
    ct=ct+1;
    end
    
    line([max(netdata{s}.decdoy) max(netdata{s}.decdoy)],[0 8],'Color','k')
 
end

figure(var+10)
subplot(3,1,2)
xlim([minday maxday])
ylim([.5 7.2])
set(gca,'Ytick',1:6,'Yticklabel',varnames)
title(sprintf('Synergistic Sources to %s',varnames{var}))
ylabel('dominant node pairs')


end

%% redundant info

for v = 1:length(plotvar)
ct=1;

var=(plotvar(v));

figure(var+10)
subplot(3,1,3)
hold on

for s = start_seg:end_seg
   
    H = entropy{s}.H_x_2(var);
 
    if H >.005
    decyear(ct) = mean(netdata{s}.decyear);
    decdoy(ct) = max(netdata{s}.decdoy);

    for i =1:6
    offset = i./50;
    StargetVar = entropy{s}.R(i,var);
    if StargetVar>1
        Wid = 12;
    elseif StargetVar>.75
        Wid = 8;
    elseif StargetVar>.25
        Wid = 5;
    else
        Wid = 1;
    end
    
    StargetVar_frac = entropy{s}.R_normbyItot(i,var);
    Wid2 = StargetVar_frac.*10;

    
    
    if StargetVar>0 %&& StargetVar==max(entropy{s}.U(:,var))
    S1lag = entropy{s}.I_dom_lag(i,var);
    StargetPair = entropy{s}.R_pair(i,var);
    S2lag = entropy{s}.I_dom_lag(StargetPair,var);
    

    figure(var+10)
    subplot(3,1,3)
    hold on
    line([decdoy(ct)-offset decdoy(ct)-offset],[i StargetPair],'LineWidth',Wid,'Color','k')
    plot(decdoy(ct)-offset, i,'.','Color',cvect(S1lag,:),'MarkerSize',Wid*5)
    plot(decdoy(ct)-offset,StargetPair,'.','Color',cvect(S2lag,:),'MarkerSize',Wid*5)

    
    end
    
    end
    
    ct=ct+1;
    end
    line([max(netdata{s}.decdoy) max(netdata{s}.decdoy)],[0 8],'Color','k')
   
 
end

figure(var+10)
subplot(3,1,3)
xlim([minday maxday])
ylim([.5 7.2])
set(gca,'Ytick',1:6,'Yticklabel',varnames)
title(sprintf('Redundant Sources to %s',varnames{var}))
ylabel('dominant node pairs')


end

