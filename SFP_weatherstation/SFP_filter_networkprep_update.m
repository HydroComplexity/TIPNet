
%SFP network prep - filter and normalize data, segment according to
%sunrise/sunset times each day

clear all
close all
clc

addpath(genpath('D:\Users\Allison\Entropy\INAS_v5_SFP'))

%load all SFP data (>800,000 data points)
load('SFP2_AllData.mat')

index = find(OrigData.Year==2015 & OrigData.data_raw_DOY>120 & OrigData.data_raw_DOY<300);
data_orig = OrigData.data_network(index,:);
decdoy = OrigData.data_raw_DOY(index);
decyear = OrigData.Year(index) + OrigData.data_raw_DOY(index)./365;
nsteps = length(index);

%leave out wind direction (index 3)
data_orig(:,3)=[];
varnames ={'RH','WS','Ta','Rg','PPT','LWet'};
nvars = 6;

%%

%replace leaf wetness index with leaf wetness raw counts, threshold = 445
lw_dat=OrigData.data_raw(index,11);
lw_dat(lw_dat<445)=445;
lw_dat(isnan(lw_dat))=445;

data_orig(:,6)=lw_dat;


figure(11)
start_DOY = 185;
end_DOY = 195;
Yr=2015;
index= find(decdoy>start_DOY & decdoy<end_DOY & floor(decyear)==Yr);

subplot(2,1,1)
plot(decdoy(index),data_orig(index,5))
ylim([0 2])

subplot(2,1,2)
plot(decdoy(index),data_orig(index,6))
ylim([400 1200])

data = data_orig;
%save('SFP_2015_DOY185_195_data.mat','data','varnames','decdoy')


%%

%identify nan values and set all nans to zero
nanflag = zeros(size(data_orig));
nanflag(isnan(data_orig))=1;
for i =1:nvars
  data_orig(nanflag==1)=0;
end

%remove PPT outliers (few extreme high values, not possible)
PPT = data_orig(:,5);
PPT(PPT>10)=0;
data_orig(:,5)=PPT;

%set RH>1 values to 1 and remove RH=0 values
RH = data_orig(:,1);
RH(RH>1)=1;
ind = find(RH==0);
for i =1:length(ind)
RH(ind(i))=RH(ind(i)-1);    
end
data_orig(:,1)=RH;

figure(1)
for i =1:nvars
    subplot(nvars,1,i)
    plot(decyear,data_orig(:,i))
    title(varnames{i})   
end

%% filter RH, Ta, and Rg (indices 1,4,5)
ind = [1 3 4];
len = [60*24/4 60*24/4 60*24/4];

data_filtered = data_orig;

ct=1;
 for i =1:3

     y = data_orig(:,ind(i));  
     [y2 a b] = ButterFiltFun('hp',y,'lambdac',len(i));

   data_filtered(:,ind(i))=y2;
   
   if ind(i)==4 %radiation, retain Rg =0 as zero
      dat = y2;
      datRg = data_orig(:,ind(i));
      dat(datRg==0)=0;
      data_filtered(:,ind(i))=dat; 
   end
   
   
     figure(2)
     subplot(length(ind),1,i)
     plot(decyear, data_filtered(:,ind(i)))
     title([varnames{ind(i)},' filtered'])
     
     figure(3)
     subplot(length(ind),1,i)
     plot(decyear, data_orig(:,ind(i))-data_filtered(:,ind(i)))
     title([varnames{ind(i)},' filtered out'])
 end
 
 

%% normalize all variables and remove outliers

rm_outliers     =[1 1 1 1 0 0];
norm            =[1 1 1 1 1 1];         

for i =1:nvars
    data_processed(:,i) = ProcessFunction(data_filtered(:,i),rm_outliers(i),norm(i));
end

figure(4)
for i =1:nvars
    subplot(nvars,1,i)
    plot(decyear,data_processed(:,i))
    title([varnames{i} ' processed'])   
end

%% segment specifically by solar radiation

Rg = data_orig(:,4);
RgIndex = zeros(size(Rg));
RgIndex(Rg~=0)=1; %daytime indicator

st=1;
dt = 30;
ct=1;

for j=1:dt:length(RgIndex)

Rgminind = RgIndex(st:st+dt-1);

if length(unique(Rgminind))>1
    
if sum(Rgminind)> (4/5)*dt %fraction of pts have Rg>0, define as day
    RgIndex(st:st+dt-1)=1;
    RgCt(ct)=1;

else            %set period as night time, set Rg values to zero
    if sum(RgIndex((st+dt-1):min((st+dt-1+100),nsteps)))< 20
    RgIndex(st:st+dt-1)=0;
    RgCt(ct)=0;
    else
    RgIndex(st:st+dt-1)=1;
    RgCt(ct)=1;
    end
end

else

RgCt(ct)=unique(Rgminind);

end

st=st+dt; 
if st+dt>length(RgIndex)
    break
end

end


segstart =1;
index = find(RgIndex==1);
segend = index(1)-1;

for seg = 1:2:10000
      
netdata{seg}.data_orig = data_orig(segstart:segend,:);
netdata{seg}.data_filtered=data_filtered(segstart:segend,:); 
netdata{seg}.data_processed=data_processed(segstart:segend,:); 
netdata{seg}.decyear = decyear(segstart:segend);
netdata{seg}.decdoy = decdoy(segstart:segend);
netdata{seg}.day_index = round(mean(RgIndex(segstart:segend)));

if netdata{seg}.day_index==0 %no radiation
    netdata{seg}.data_processed(:,5)=0;
end

netdata{seg}.data_process2 = netdata{seg}.data_processed;
seglen(seg)=length(segstart:segend);
netdata{seg}.seglength = seglen(seg);
yr(seg) = mean(decyear(segstart:segend));

segstart = segend+1; %beginning of next segment

%index of next value where Rg changes
RgshortIndex = RgIndex(segstart:end); %remove finished values
index = find(RgshortIndex~=RgIndex(segstart)); %index(1) is where Rg shifts again
if length(index)>1
segend = segstart+floor((index(1)/2));
segend2 = segstart+index(1)-1;
else
segend = length(RgIndex);
end

if segend >= nsteps || segstart>=nsteps
    break
end

netdata{seg+1}.data_orig = data_orig(segstart:segend,:);
netdata{seg+1}.data_filtered=data_filtered(segstart:segend,:); 
netdata{seg+1}.data_processed=data_processed(segstart:segend,:); 
netdata{seg+1}.decyear = decyear(segstart:segend);
netdata{seg+1}.decdoy = decdoy(segstart:segend);
netdata{seg+1}.day_index = round(mean(RgIndex(segstart:segend)));

if netdata{seg+1}.day_index==0 %no radiation
    netdata{seg+1}.data_processed(:,5)=0;
end

netdata{seg+1}.data_process2 = netdata{seg+1}.data_processed;
seglen(seg+1)=length(segstart:segend);
yr(seg+1) = mean(decyear(segstart:segend));
netdata{seg+1}.seglength = seglen(seg+1);

segstart=segend+1;
segend=segend2;

if segend >= nsteps || segstart>=nsteps
    break
end 
    
end

figure(6)
plot(yr,seglen)

filename = 'SFPproject_032717.mat';
save(filename,'netdata','OrigData','varnames','data_processed')

%% set up mi structure to run network code directly

mi.nvars = 6;
mi.nSteps =size(OrigData.data_raw,1);
mi.segsteps = 400;
mi.nsegs =  length(netdata);
mi.netopt = 3;
mi.nTests = 100;
mi.N = 35;
mi.lagvect = [1:10 15:5:60];
mi.nlags = length(mi.lagvect);

mi.bin_scheme = 'global';
mi.maxS = 20;
mi.Range = [min(data_processed); max(data_processed)];
mi.method = 'KDE';
mi.parallel_opt =1;
mi.ZeroLagOpt = 0;
mi.NoSelfOpt = 0;
mi.DomNormOpt = 0;

for i =1:mi.nvars
    
  dat = data_processed(:,i);
  R = range(dat);
  n = mi.segsteps;
  if i>5
      n=50; %more zero values, fewer actual values for PPT and LW
  end
  
  dat_nonzero = dat;
  dat_nonzero(dat==0)=[];
  vardat = var(dat_nonzero);
  
  mi.KDEparams.h1_seg_orig(i)= 1.06 * n^(-1/5) .* vardat./R;
  mi.KDEparams.h2_seg_orig(i)= 1.77 * n^(-1/6) .* vardat./R;
  mi.KDEparams.h3_seg_orig(i) = 2.78 * n^(-1/7) .* vardat./R;
end

mi.KDEparams.h1 = mi.KDEparams.h1_seg_orig;
mi.KDEparams.h2 = mi.KDEparams.h2_seg_orig;
mi.KDEparams.h3 = mi.KDEparams.h3_seg_orig;

save(filename,'mi','-append')

%% run entropy code
%clear all
%close all
% clc
% 

parfor_progress(mi.nsegs);
for i =1:mi.nsegs
    entropy{i}=1;
end

parfor i=1:mi.nsegs
    sprintf('%d\n',i)
    dat = netdata{i}.data_process2; 
    entropy{i}=EntropyFun_SFPUpdate(mi,dat,i);
    parfor_progress;
   sprintf('%d\n',i)
end

parfor_progress(0);

save(filename,'entropy','mi','netdata','OrigData','-append');

%AllStats = TotalStatsFunctionGUI(netdata,mi,entropy);
%save(filename,'AllStats','-append')