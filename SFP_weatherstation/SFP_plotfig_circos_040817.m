
%SFP network plotting of frequent pairs of S and R sources

clear all
close all
clc

addpath(genpath('D:\Users\Allison\Entropy\INAS_v5_SFP'))

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

close all

varvect = 1:6;
nvars = length(varvect);

dayind =[0 0 0 1 1 1]; %night then day
wetind =[1 0 2 1 0 2]; %1=dew wet, 2=rain, 3 = dry
numsegs=zeros(1,6);

I_lag = zeros(nvars,nvars,7);
I_lagall = zeros(nvars,nvars);


%time scale of interactions

ct=ones(nvars,nvars,7);
ctall = ones(nvars,nvars);
numdd = ones(1,7);

for s=start_seg:end_seg
    
    SumPPT = sum(netdata{s}.data_orig(:,5)>0);
    SumLWet = sum(netdata{s}.data_processed(:,6)>0);
    
    if  SumPPT>5 && SumLWet>100 && netdata{s}.day_index==1
        dd = 1; %wet
        numdd(dd)=numdd(dd)+1;
    elseif SumLWet==0 && SumPPT==0 && netdata{s}.day_index==1
        dd=2; %day dry
        numdd(dd)=numdd(dd)+1;
    elseif SumLWet>100 && SumPPT==0 && netdata{s}.day_index==1
        dd=3;  %day dew
        numdd(dd)=numdd(dd)+1;
    elseif SumPPT>5 && SumLWet>100 && netdata{s}.day_index==0
        dd = 4; %wet
        numdd(dd)=numdd(dd)+1;
    elseif SumLWet==0 && SumPPT==0 && netdata{s}.day_index==0
        dd=5; %
        numdd(dd)=numdd(dd)+1;
    elseif SumLWet>100 && SumPPT==0 && netdata{s}.day_index==0
        dd=6;
        numdd(dd)=numdd(dd)+1;
    else
        dd=7;
    end

    for v1=1:length(varvect)
        for v2=1:length(varvect)
            if ~isnan(entropy{s}.I_dom_lag(v1,v2))
               
                I_lag(v1,v2,dd)=I_lag(v1,v2,dd)+entropy{s}.I_dom_lag(v1,v2);
                ct(v1,v2,dd)=ct(v1,v2,dd)+1;

                I_lagall(v1,v2)=I_lagall(v1,v2)+entropy{s}.I_dom_lag(v1,v2);
                ctall(v1,v2)=ctall(v1,v2)+1;
            end
        end
    end
    % I_lag(v1,v2)=I_lag(v1,v2)./ct;
end


%%

for dd=1:6 %night/day index
pair =1;
numpairs = nchoosek(nvars-1,2);
pairvars = zeros(numpairs,2);
RS_count = zeros(numpairs,nvars*3);
RS_strength = zeros(numpairs,nvars*3);
U_count = zeros(nvars*3,nvars*3);
U_strength = zeros(nvars*3,nvars*3);

for v1 = 1:length(varvect)
    n1 = nvars-varvect(v1)+1;
    for v2 = 1:(n1-1)
        n2 = varvect(v2);
        
        if n1==5 || n2==5 %skip PPT
            continue
        end

        pairvars(pair,:)=[n1 n2];
        pairstring{pair} = sprintf('%s%s',varnames{n1},varnames{n2});
        
        tarct=1;
        for tar = 1:nvars
            
            if tar == 4 || tar ==5 %skip Rg as a target
                tarct=tarct+1;
                continue
            end
                      
            for s =start_seg:end_seg
                
                SumPPT = sum(netdata{s}.data_orig(:,5)>0);
                SumLWet = sum(netdata{s}.data_processed(:,6)>0);
                
                if netdata{s}.day_index ~=dayind(dd) %separate night and day times
                    continue
                elseif (wetind(dd)==1 && SumLWet<100) || (wetind(dd)==1 && SumPPT>0)
                    continue
                elseif wetind(dd)==2 && SumPPT<5
                    continue
                elseif (wetind(dd)==0 && SumLWet>0) || (wetind(dd)==0 && SumPPT>0)
                    continue
                elseif dayind(dd)==0 && sum(entropy{s}.R(4,:))>0
                    continue
                end
                
                if v1==1 &&v2==2 &&tar==1
                    numsegs(dd)=numsegs(dd)+1;             
                end
      
                S1 = entropy{s}.S(n1,tar);
                S2 = entropy{s}.S(n2,tar);
                R1 = entropy{s}.R(n1,tar);
                R2 = entropy{s}.R(n2,tar);
                Itot1 = entropy{s}.Itot(n1,tar);
                Itot2 = entropy{s}.Itot(n2,tar);
                   
                if R1>0 && entropy{s}.R_pair(n1,tar)==n2 %&& R1 == max(entropy{s}.R(:,tar))
                    RS_count(pair,tarct) = RS_count(pair,tarct)+1;
                    RS_strength(pair,tarct) = RS_strength(pair,tarct)+entropy{s}.R(n1,tar);
                    
                    U_count(n1+2*nvars,tarct+2*nvars) = U_count(n1+2*nvars,tarct+2*nvars)+1;
                    U_strength(n1+2*nvars,tarct+2*nvars) = U_strength(n1+2*nvars,tarct+2*nvars)+entropy{s}.U(n1,tar); 
                    
                elseif R2>0 && entropy{s}.R_pair(n2,tar)==n1 %&& R2 == max(entropy{s}.R(:,tar))
                    RS_count(pair,tarct) = RS_count(pair,tarct)+1;
                    RS_strength(pair,tarct) = RS_strength(pair,tarct)+entropy{s}.R(n2,tar);  
                    
                    U_count(n2+2*nvars,tarct+2*nvars) = U_count(n2+2*nvars,tarct+2*nvars)+1;
                    U_strength(n2+2*nvars,tarct+2*nvars) = U_strength(n2+2*nvars,tarct+2*nvars)+entropy{s}.U(n2,tar);
                    
                end
                
                if S1>0 && entropy{s}.S_pair(n1,tar)==n2 %&& S1 == max(entropy{s}.S(:,tar))
                    RS_count(pair,tarct+nvars) = RS_count(pair,tarct+nvars)+1;
                    RS_strength(pair,tarct+nvars) = RS_strength(pair,tarct+nvars)+entropy{s}.S(n1,tar);     
                elseif S2>0 && entropy{s}.S_pair(n2,tar)==n1 %&& S2 == max(entropy{s}.S(:,tar))
                    RS_count(pair,tarct+nvars) = RS_count(pair,tarct+nvars)+1;
                    RS_strength(pair,tarct+nvars) = RS_strength(pair,tarct+nvars)+entropy{s}.S(n2,tar);
                end
                
            end
         
        tarct = tarct+1;
        end
        pair = pair+1;  
    end
end


RS_percent = RS_count./numsegs(dd);
RS_strength = RS_strength./RS_count;
RS_strength(RS_count==0)=0;

U_percent = U_count./numsegs(dd);
U_strength = U_strength./U_count;
U_strength(U_count==0)=0;

RS_percent(RS_percent<.25)=0;
%RS_percent = RS_percent.*RS_strength;
RS_strength(RS_percent==0)=0;

U_percent(U_percent<.25)=0;
%U_percent = U_percent.*U_strength;
U_strength(U_percent==0)=0;

if dd==1 %night time networks
    RS_strength_night1_c = zeros(nvars*3+numpairs);
    RS_strength_night1_c((nvars*3+1):end,(numpairs+1):end)=RS_strength;
    RS_strength_night1_c(1:nvars*3,(numpairs+1):end)=U_strength;
    
    RS_percent_night1_c = zeros(nvars*3+numpairs);
    RS_percent_night1_c((nvars*3+1):end,(numpairs+1):end)=RS_percent; 
    RS_percent_night1_c(1:nvars*3,(numpairs+1):end)=U_percent;
    
elseif dd==2
    RS_strength_night2_c = zeros(nvars*3+numpairs);
    RS_strength_night2_c((nvars*3+1):end,(numpairs+1):end)=RS_strength;
    RS_strength_night2_c(1:nvars*3,(numpairs+1):end)=U_strength;
    
    RS_percent_night2_c = zeros(nvars*3+numpairs);
    RS_percent_night2_c((nvars*3+1):end,(numpairs+1):end)=RS_percent;
    RS_percent_night2_c(1:nvars*3,(numpairs+1):end)=U_percent;
elseif dd==3
    RS_strength_night3_c = zeros(nvars*3+numpairs);
    RS_strength_night3_c((nvars*3+1):end,(numpairs+1):end)=RS_strength;
    RS_strength_night3_c(1:nvars*3,(numpairs+1):end)=U_strength;
    
    RS_percent_night3_c = zeros(nvars*3+numpairs);
    RS_percent_night3_c((nvars*3+1):end,(numpairs+1):end)=RS_percent;
    RS_percent_night3_c(1:nvars*3,(numpairs+1):end)=U_percent;
elseif dd==4
    RS_strength_day1_c = zeros(nvars*3+numpairs);
    RS_strength_day1_c((nvars*3+1):end,(numpairs+1):end)=RS_strength;
    RS_strength_day1_c(1:nvars*3,(numpairs+1):end)=U_strength;
    
    RS_percent_day1_c = zeros(nvars*3+numpairs);
    RS_percent_day1_c((nvars*3+1):end,(numpairs+1):end)=RS_percent;
    RS_percent_day1_c(1:nvars*3,(numpairs+1):end)=U_percent;
elseif dd==5
    RS_strength_day2_c = zeros(nvars*3+numpairs);
    RS_strength_day2_c((nvars*3+1):end,(numpairs+1):end)=RS_strength;
    RS_strength_day2_c(1:nvars*3,(numpairs+1):end)=U_strength;
    
    RS_percent_day2_c = zeros(nvars*3+numpairs);
    RS_percent_day2_c((nvars*3+1):end,(numpairs+1):end)=RS_percent;
    RS_percent_day2_c(1:nvars*3,(numpairs+1):end)=U_percent;
elseif dd==6
    RS_strength_day3_c = zeros(nvars*3+numpairs);
    RS_strength_day3_c((nvars*3+1):end,(numpairs+1):end)=RS_strength;
    RS_strength_day3_c(1:nvars*3,(numpairs+1):end)=U_strength;
    
    RS_percent_day3_c = zeros(nvars*3+numpairs);
    RS_percent_day3_c((nvars*3+1):end,(numpairs+1):end)=RS_percent;
    RS_percent_day3_c(1:nvars*3,(numpairs+1):end)=U_percent;
end


end

nvars = nvars*3;
varnames={'RH_R','WS_R','Ta_R','Rg_R','PPT_R','LWet_R','RH_S','WS_S','Ta_S','Rg_S',...
    'PPT_S','LWet_S','RH_U','WS_U','Ta_U'...
    ,'Rg_U','PPT_U','LWet_U'};


%% create text files for CIRCOS plots

Matrix = RS_percent_day1_c;

foldername = 'TextFiles';
addpath(foldername);
if exist(foldername)~=7
    mkdir(foldername);
end
for ff=1:6
    
if ff==1
 fileID = fopen([foldername,'/RS_strength_day1.txt'],'w+');  
 Matrix = RS_percent_day1_c;
elseif ff==2
 fileID = fopen([foldername,'/RS_strength_day2.txt'],'w+');
 Matrix = RS_percent_day2_c;
elseif ff==3
  fileID = fopen([foldername,'/RS_strength_day3.txt'],'w+');
  Matrix = RS_percent_day3_c;
elseif ff==4
   fileID = fopen([foldername,'/RS_strength_night1.txt'],'w+');
   Matrix = RS_percent_night1_c;
elseif ff==5
  fileID =  fopen([foldername,'/RS_strength_night2.txt'],'w+');
  Matrix = RS_percent_night2_c;
else
  fileID =  fopen([foldername,'/RS_strength_night3.txt'],'w+');
  Matrix = RS_percent_night3_c;
end

fprintf(fileID,'label label ');

for i =1:nvars+numpairs
fprintf(fileID,'%d ',i);
end

fprintf(fileID,'\n');




fprintf(fileID,'label label ');


for i =1:numpairs   
    fprintf(fileID,'%s ',pairstring{i});
end

for i =1:nvars
    fprintf(fileID,'%s ',varnames{i});

end

fprintf(fileID,'\n');

col_order = [1 7 13 2 8 14 3 9 15 4 10 16 5 11 17 6 12 18:(nvars+numpairs)]; 

ci=200; %gray color
col = [255 0 0; 255 120 0; 255 190 0;255 255 0;  155 255 47; 34 139 34;...
    0 245 255;0 0 255;127 0 255;153 50 204]; %vector of colors for variable pairs

 for i = 1:nvars+numpairs
     if i <=6
         fprintf(fileID,'%d,%d,%d %s ',0 ,0,0, varnames{i});
     elseif i <=12
         fprintf(fileID,'%d,%d,%d %s ',255, 255, 255, varnames{i});
     elseif i<=18
         fprintf(fileID,'%d,%d,%d %s ',100, 100, 100, varnames{i});
     else
         fprintf(fileID,'%d,%d,%d %s ',col(i-nvars,1),col(i-nvars,2),col(i-nvars,3), pairstring{i-nvars});
     end
     
     for j = 1:nvars+numpairs
         fprintf(fileID,'%d  ',round(Matrix(i,j)*100));
     end
     fprintf(fileID,'\n');

 end
end
