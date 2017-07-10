
%plots from Goodwell and Kumar 2017 (TIPNets)
%plots of dominant average time scales (minutes) and differences for
%weather categories, bar plot of total information and proportions of U, R,
%and S for different weather categories

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

%%


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

I_lagvar1 = zeros(nvars,nvars,5);
I_lagvar2 = zeros(nvars,nvars,5);
I_lagvar3 = zeros(nvars,nvars,5);

%time scale of interactions

ct=zeros(nvars,nvars,7);
ctall = zeros(nvars,nvars);
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
    
    day_ind(s)=dd;
    % I_lag(v1,v2)=I_lag(v1,v2)./ct;
end

I_lagall=I_lagall./ctall;
I_lag = I_lag./ct;

I_lagdiff1 = round(I_lag(1:4,1:3,1)-I_lagall(1:4,1:3));
I_lagdiff2 = round(I_lag(1:4,1:3,2)-I_lagall(1:4,1:3));
I_lagdiff3 = round(I_lag(1:4,1:3,3)-I_lagall(1:4,1:3));
I_lagdiff4 = round(I_lag(1:3,1:3,4)-I_lagall(1:3,1:3));
I_lagdiff5 = round(I_lag(1:3,1:3,5)-I_lagall(1:3,1:3));
I_lagdiff6 = round(I_lag(1:3,1:3,6)-I_lagall(1:3,1:3));

%check if differences in timescale should be omitted (less than 10 in
%category)
minct = 20;
omit1 = reshape(ct(:,:,1),nvars,nvars);
omit1(omit1<minct)=1;
omit1(omit1>=minct)=0;
omit2 = reshape(ct(:,:,2),nvars,nvars);
omit2(omit2<minct)=1;
omit2(omit2>=minct)=0;
omit3 = reshape(ct(:,:,3),nvars,nvars);
omit3(omit3<minct)=1;
omit3(omit3>=minct)=0;
omit4 = reshape(ct(:,:,4),nvars,nvars);
omit4(omit4<minct)=1;
omit4(omit4>=minct)=0;
omit5 = reshape(ct(:,:,5),nvars,nvars);
omit5(omit5<minct)=1;
omit5(omit5>=minct)=0;
omit6 = reshape(ct(:,:,6),nvars,nvars);
omit6(omit6<minct)=1;
omit6(omit6>=minct)=0;

bluewhitered = [.4 0 0; .8 0 0; 1 .3 .3; 1 .6 .6; 1 1 1; .6 .6 1; .3 .3 1; 0 0 .8; 0 0 .4];

figure(4)
h1 = axes('units','centimeters','position',[1 1 5 3]);
imagesc(I_lagall);
colormap(flipud(jet(7)))
set(h1,'Ytick',1:6,'Yticklabel',varnames,'Xtick',1:6,'Xticklabel',varnames)
colorbar
set(h1,'FontSize',9)
caxis([0 35])
title('Dominant Lag times')
colorbar

print('-depsc','-painters','-loose','Fig_TimeDayAvg')

figure(5)
h1 = axes('units','centimeters','position',[1 5 4 3]);
imagesc(I_lagdiff1);
set(h1,'Ytick',1:6,'Yticklabel',varnames,'Xtick',1:6,'Xticklabel',varnames)
colormap(bluewhitered)
colorbar
set(h1,'FontSize',9)
caxis([-13.5 13.5])
title('Rainy Day - Avg')

h1 = axes('units','centimeters','position',[6 5 4 3]);
imagesc(I_lagdiff2);
colormap(bluewhitered)
colorbar
set(h1,'Ytick',1:6,'Yticklabel',varnames,'Xtick',1:6,'Xticklabel',varnames)
colorbar
set(h1,'FontSize',9)
caxis([-13.5 13.5])
title('Dry Day - Avg')

h1 = axes('units','centimeters','position',[11 5 4 3]);
imagesc(I_lagdiff3);
colormap(bluewhitered)
colorbar
set(h1,'Ytick',1:6,'Yticklabel',varnames,'Xtick',1:6,'Xticklabel',varnames)
colorbar
set(h1,'FontSize',9)
caxis([-13.5 13.5])
title('Dew Day - Avg')

h1 = axes('units','centimeters','position',[1 1 4 3]);
imagesc(I_lagdiff4);
set(h1,'Ytick',1:6,'Yticklabel',varnames,'Xtick',1:6,'Xticklabel',varnames)
colormap(bluewhitered)
colorbar
set(h1,'FontSize',9)
caxis([-13.5 13.5])
title('Rainy Night - Avg')

h1 = axes('units','centimeters','position',[6 1 4 3]);
imagesc(I_lagdiff5);
colormap(bluewhitered)
colorbar
set(h1,'Ytick',1:6,'Yticklabel',varnames,'Xtick',1:6,'Xticklabel',varnames)
colorbar
set(h1,'FontSize',9)
caxis([-13.5 13.5])
title('Dry Night - Avg')

h1 = axes('units','centimeters','position',[11 1 4 3]);
imagesc(I_lagdiff6);
colormap(bluewhitered)
colorbar
set(h1,'Ytick',1:6,'Yticklabel',varnames,'Xtick',1:6,'Xticklabel',varnames)
colorbar
set(h1,'FontSize',9)
caxis([-13.5 13.5])
title('Dew Night - Avg')

print('-depsc','-painters','-loose','Fig_TimeDiff')

%%

for dd=1:6 %night/day index
pair =1;
numpairs = nchoosek(nvars,2);
pairvars = zeros(numpairs,2);
RS_count = zeros(numpairs,nvars*3);
RS_strength = zeros(numpairs,nvars*3);
U_count = zeros(nvars*3,nvars*3);
U_strength = zeros(nvars*3,nvars*3);


for v1 = 1:length(varvect)
    n1 = nvars-varvect(v1)+1;
    for v2 = 1:(n1-1)
        n2 = varvect(v2);
        
        if n1==5 % || n2==5 %skip PPT
            continue
        end

        pairvars(pair,:)=[n1 n2];
        pairstring{pair} = sprintf('%s%s',varnames{n1},varnames{n2});
        
        tarct=1;
        for tar = 1:nvars
            
            if tar == 4 %|| tar ==5 %skip Rg as a target
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


totalinfo(dd)=sum(sum(RS_strength))+sum(sum(U_strength));
totalu(dd) = sum(sum(U_strength));
totalr(dd) = sum(sum(RS_strength(:,1:nvars)));
totals(dd) = sum(sum(RS_strength(:,(nvars+1):end)));

end

nvars = nvars*3;


%% plot bar chart of relative information measures

U_frac = totalu./totalinfo;
S_frac = totals./totalinfo;
R_frac = totalr./totalinfo;

dat = [U_frac; S_frac; R_frac];
dat2 = [totalu; totals; totalr];

figure(11)
subplot(2,1,1)
hh = bar(dat2','stacked');


set(gca,'Xticklabel',{'wet night','dry night','rain night','wet day','dry day','rain day'})
ylabel('bits')
title('Total avg Contributions')
P=findobj(gca,'type','patch');
C=[0 0 0; 1 1 1; .5 .5 .5]; % make a colors list 
for n=1:length(P) 
set(P(n),'facecolor',C(n,:));
end
legend('U','S','R')

subplot(2,1,2)
hh=bar(dat','stacked');
set(gca,'Xticklabel',{'wet night','dry night','rain night','wet day','dry day','rain day'})
legend(gca,'U','S','R')
title('Relative Contributions')
ylim([0 1])
ylabel('fraction')
P=findobj(gca,'type','patch');
C=[0 0 0; 1 1 1; .5 .5 .5]; % make a colors list 
for n=1:length(P) 
set(P(n),'facecolor',C(n,:));
end


