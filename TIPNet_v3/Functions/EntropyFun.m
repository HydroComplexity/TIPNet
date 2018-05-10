function [entropy] = EntropyFun(mi,Data,seg)
%Function to compute various entropy, transfer entropy, information decomposition measures
%Allison Goodwell, June 2014
%August 2014: also calculate lagged mutual information
%Jan 2015: updated method of KDE
%Jan 2015: alter method of transfer entropy: compute T based on
%conditioning lag of min[H(Yt)|H(Yt-lag))]
%Feb 2015: bin_scheme: 'local' or 'global', Range: min and max data values
%Feb 2015: changed kd combined function to normalize variables
%June 2015: changed KDE and info measures code (to also compute
%redundant/synergistic information
%July 2015: save certain computed pdfs according to [trans, rec, lag]
%vector
%July 2015: consider multiple sources as contributors to target
%maxS = number of sources to any single target varible to compute unique
%information
%August 2015: modify code to reduce number of computations to redundancy:
%first eliminate redundant links from each source node (with multiple
%signifcant lags - many lags likely to be completely redundant)
%September 2015: added cumulative synergy, changed outputs and inputs to
%structures: mi has all parameters, data is matrix of data, entropy
%is structure with entropy results and pdfs according to what should be
%saved
%September 8, 2015: implement alternate version of redundancy (also look
%for correlation between sources)
%September 30, 2015: updated compute_info_measures for output as structure
%January 18, 2016: updates for GUI - input h (smoothing parameter) to KDE
%computation, option for method in mi - fixed binning or KDE
%February 2016: update options: mi.netopt (1=full network, 2 = H and I
%only, 3 = H only (only entropy of nodes)
%February 2016: alter redundancy: 2 redundancy matrices: R_tau is
%redundancy between a single source to a target node at different time
%lags, R_T is redundancy between dominant link from each source to a target

%Update: altered version to consider zero-lag I as possible dominant link
%for T/I, other measures, depending on mi.ZeroLagOpt

%Update 6/10/16: mi.DomNormOpt = 0 for classifying dominant link as
%non-normalized (default), = 1 for normalized (previous method)
%changed several variable names

nvars = mi.nvars;
nTests = mi.nTests;
N = mi.N;
lagvect = mi.lagvect;
bin_scheme=mi.bin_scheme;
Range = mi.Range;
method = mi.method;
opt = mi.netopt;
z_opt = mi.ZeroLagOpt;
DomNormOpt = 0;
z_effect = mi.DataPrep.Z_effect;


ndata = size(Data,1);

if z_opt==1 %include zero in lagvect
lagvect =[0 lagvect];
end

for i =1:nvars
  
    dat = Data(:,i);
    Ra = range(dat);
    varterm = var(dat);
    varterm(isnan(varterm))=0;
    
h1D(i)= 1.06 * ndata^(-1/5) .* varterm./Ra; 
h2D(i)= 1.77 * ndata^(-1/6) .* varterm./Ra; 
h3D(i) = 2.78 * ndata^(-1/7) .* varterm./Ra; 
end

nlags = length(lagvect);
lim = eps*10; % machine limit of matlab 
N_sources = zeros(nvars,1);
N_sourcenodes = zeros(nvars,1);
TE = zeros(nlags,nvars,nvars);
TE_normbyItot = TE;
I_lags = zeros(nlags,nvars,nvars);
I_normbyH = zeros(nlags,nvars,nvars);
I_inst = zeros(nvars,nvars);
I_inst_normbyH = zeros(nvars,nvars);
I_inst_sig = zeros(nvars,nvars);
H_x_1 = zeros(nvars,1);
H_x_2 = zeros(nvars,1);

TE_T_normbyItot = zeros(nvars,nvars);

I_pvalue = zeros(nvars);
Itot = zeros(nvars);

S_allpairs = zeros(nchoosek(nvars,2),nvars);
R_allpairs = zeros(nchoosek(nvars,2),nvars);
Itot_allpairs = zeros(nchoosek(nvars,2),nvars);
TI_allpairs = zeros(nchoosek(nvars,2),nvars);
I_ind1_allpairs = zeros(nchoosek(nvars,2),nvars);
I_ind2_allpairs = zeros(nchoosek(nvars,2),nvars);

S = zeros(nvars);
R = S;
U = S;
S_pair = S;  %synergy buddies
R_pair = S;  %redundant buddies
TI_pair = S; %T/I buddies (second source which minimizes T/I for given source)
Itot_pair = S;


TI = zeros(nvars);
TE_T = TI;
Itot_T = zeros(nvars);
I_dom = zeros(nvars);
I_dom_normbyH = zeros(nvars);
I_tau = nan(nvars);

%% flag variables with constant values or NaN values, don't use
flag = zeros(1,nvars);
for Source=1:nvars
    X= Data(:,Source);
    numNaNs = sum(isnan(X));
    num9999 = sum(X==-9999);
    %check for all zero or all constant values
    if std(X)<lim || numNaNs >0 || num9999 >0
        flag(Source)=1;
    end
end

%% compute 1D entropy values for each node
for Source = 1:nvars
    
    if flag(Source)==1
        continue
    end
    X = Data(:,Source);
    
    pdf = compute_pdfGUI(X,N,bin_scheme, Range(:,Source),method,z_effect(Source+1));
    info = compute_info_measures(pdf);
    H_x_1(Source) = info.Hx;
end

if opt==1                                                   %done, exit function
    entropy.H_x_1 = H_x_1;
    return
end

%% compute MI and lagged MI for all lags and pairs

timer=timebar(sprintf('MI for segment %d',seg),'Network Timer');
ct=0;
timebar(timer,ct/nvars^2);

for Source = 1:nvars %transmitters
    
    X = Data(:,Source);
    
    if flag(Source)==1
        ct=ct+nvars;
        continue
    end
    
    for Target = 1:nvars %receivers
        %fprintf('trans = %d, rec = %d\n',Source,Target)
        
        Y = Data(:,Target);
        if flag(Target)==1
            ct=ct+1;
            continue
        end
   
        for t=1:nlags %may or may not include zero lag           

            lagt = lagvect(t);

            nTuples = ndata-lagt-1;
            tar_start=1+lagt;
            tarlag_start = max(1,lagt);
            
            Svar=X(1:nTuples);        %Leading Node S (lag tau earlier than present)
            Tar=Y(tar_start:tar_start+nTuples-1);        %Led Node Tar (one timestep in future)
            TarLag=Y(tarlag_start:tarlag_start+nTuples-1);        %current node TarLag
            
            %first compute mutual information for all lags
            % fprintf('computing lagged I Source %d Target %d lag %d\n',Source,Target,lagt)
            if std(Svar)>lim && std(Tar)>lim
                pdf = compute_pdfGUI([Svar Tar],N,bin_scheme, [Range(:,Source)  Range(:,Target)],method,[z_effect(Source+1) z_effect(Target+1)]);
                info = compute_info_measures(pdf);
                I_lags(t,Source,Target)=info.I;
                Htemp = min(info.Hx1,info.Hx2);
                
            %shuffle to compute significance
            if nTests>0
                if t==1 %do shuffled significance testing for zero lag only
                    I_shuff=zeros(1,nTests);
                    for test=1:nTests
                        Tar_shuff = randsample(Tar,length(Tar));
                        pdfshuff = compute_pdfGUI([Svar Tar_shuff],N,bin_scheme, ...
                            [Range(:,Source)  Range(:,Target)],method,[z_effect(Source+1) z_effect(Target+1)]);
                        infoshuff = compute_info_measures(pdfshuff);
                        I_shuff(test) = infoshuff.I;
                    end
                    I_shuff_sig = mean(I_shuff)+4*std(I_shuff);
                    [h,p,ci,stats] = ttest((info.I-I_shuff)./std(I_shuff));
                end

                I_pvalue(Source,Target)=p;
                I_inst_sig(Source,Target)=I_shuff_sig;
 
                if info.I-I_shuff_sig< 0
                I_lags(t,Source,Target) = 0;
                end    
            end
            
            I_normbyH(t,Source,Target)=I_lags(t,Source,Target)./Htemp;
            
            % for significant lagged I, also compute Transfer Entropy TE
            if I_lags(t,Source,Target)>0 && opt==3 
                pdf = compute_pdfGUI([Svar TarLag Tar],N,bin_scheme,...
                    [Range(:,Source) Range(:,Target) Range(:,Target)],method,[z_effect(Source+1) z_effect(Target+1)]);
                info = compute_info_measures(pdf);
                TE(t,Source,Target)=info.T;
                if info.Itot>0
                TE_normbyItot(t,Source,Target) = info.T/info.Itot; %normalized by total information shared
                end
               
            end
            
            end
            
        end                                                         %lags
        
        if z_opt ==1 %set I_inst to first index of I_lags
            I_inst(Source,Target) = I_lags(1,Source,Target);
            I_inst_normbyH(Source,Target) = I_normbyH(1,Source,Target);
        else   %otherwise compute instantaneous mutual info separately
            pdf = compute_pdfGUI([X Y],N,bin_scheme, [Range(:,Source)  Range(:,Target)],method, [z_effect(Source+1) z_effect(Target+1)]);
            info = compute_info_measures(pdf);
            I_inst(Source,Target)=info.I;
            if nTests>0
                I_shuff=zeros(1,nTests);
                for test=1:nTests
                    Tar_shuff = randsample(Y,length(Y));
                    
                    pdfshuff = compute_pdfGUI([X Tar_shuff],N,bin_scheme, ...
                        [Range(:,Source)  Range(:,Target)],method, [z_effect(Source+1) z_effect(Target+1)]);
                    infoshuff = compute_info_measures(pdfshuff);
                    I_shuff(test) = infoshuff.I;
                end
                I_shuff_sig = mean(I_shuff)+4*std(I_shuff);
                
                I_inst_sig(Source,Target)=I_shuff_sig;
               
                if I_shuff_sig - info.I > 0
                
                I_inst(Source,Target) =0; 
                
                else
                I_inst_normbyH(Source,Target) = I_inst(Source,Target)./min(info.Hx1,info.Hx2);
                end
                
                
            end
            
        end
        
        ct=ct+1;
        timebar(timer,ct/nvars^2)
    end                                                             %receivers
end                                                                 %transmitters

%re-define H_x as I(X;X) (from the 2D pdf instead of 1D)
for n =1:nvars
H_x_2(n) = I_inst(n,n);
end

I_normbyH(isnan(I_normbyH))=0;
I_inst_normbyH(isnan(I_inst_normbyH))=0;

%matrix of dominant links
for Target =1:nvars
    for Source = 1:nvars
        
        vect_norm = reshape(I_normbyH(:,Source,Target),1,nlags); %normalized I (bits/bit)
        vect_nonnorm = reshape(I_lags(:,Source,Target),1,nlags);   %non-normalized I (bits)
        
        if Source==Target && mi.ZeroLagOpt ==1
            vect_norm(1)=0; %don't allow for zero-lag self links
            vect_nonnorm(1)=0;
        end
        
        if Source==Target && mi.NoSelfOpt==1 %skip self links if NoSelfOpt=1
            continue
        end

        if max(vect_norm)>0 && DomNormOpt == 1 %define dominant links based on normalized values
            [I_dom_normbyH(Source,Target), ind_lag] = max(vect_norm);    %maximum strength
            I_tau(Source,Target) = lagvect(ind_lag(1));
            TE_T(Source,Target) = TE(ind_lag(1),Source,Target);
            TE_T_normbyItot(Source,Target) = TE_normbyItot(ind_lag(1),Source,Target);

            
            I_dom(Source,Target) = I_lags(ind_lag,Source,Target);
        elseif max(vect_nonnorm>0) %define dominant links based on I in bits (non-normalized)
            [I_dom(Source,Target), ind_lag] = max(vect_nonnorm);    %maximum strength
            I_tau(Source,Target) = lagvect(ind_lag(1));
            TE_T(Source,Target) = TE(ind_lag(1),Source,Target);
            TE_T_normbyItot(Source,Target) = TE_normbyItot(ind_lag(1),Source,Target);
            
            I_dom_normbyH(Source,Target) = I_normbyH(ind_lag,Source,Target);
        end            
    end   
end

if DomNormOpt==1
    IDOM = I_dom_normbyH;
else
    IDOM = I_dom;
end

close(timer)

if opt==2                                                   %done, exit function
    entropy.H_x_1 =             H_x_1;
    entropy.H_x_2 =             H_x_2;
    entropy.I_tau =             I_lags;
    entropy.I_normbyH =         I_normbyH;
    entropy.I_inst =            I_inst;
    entropy.I_inst_normbyH =    I_inst_normbyH;
    entropy.I_dom =             I_dom;
    entropy.I_dom_normbyH =     I_dom_normbyH;
    entropy.I_dom_lag =         I_tau;
    return
end

%% use lagged MI to compute synergies/redundancies between sources

timer=timebar(sprintf('Redundancy for segment %d',seg),'Network Timer');
timebar(timer,0/nvars);

%if option to not consider self links on - omit all self links
if mi.NoSelfOpt ==1
    for i =1:nvars
        I_lags(:,i,i)=0;
    end
end

pairct = 1; %count of pairs

%pairwise analysis - determine S and R to each target from each source pair
for v1 = 1:nvars
    n1 = nvars-v1+1;
    for v2 = 1:(n1-1)
        n2 = v2;
        pairvars(pairct,:)=[n1 n2];
        %pairstring{pairct} = sprintf('%s%s',varnames{n1},varnames{n2});
        
        for Target = 1:nvars
            if flag(Target)==1 || sum(sum(I_lags(:,:,Target)))==0 ||...
                    flag(n1)==1 || flag(n2)==1
                continue        %go to next target node if Target has no sources
            elseif I_dom(n1,Target)==0 || I_dom(n2,Target)==0
                continue        %go to next target if one source does not provide info
            end
            
            %both sources provide info to target: compute total I and partition
            lag1 = I_tau(n1,Target);
            lag2 = I_tau(n2,Target);
            nTuples = ndata-max(lag1,lag2)-1;
            
            if lag1>lag2
                s1_start=1;
                s2_start=1+lag1-lag2;
                tar_start=1+lag1;
            elseif lag2>lag1
                s1_start=1+lag2-lag1;
                s2_start=1;
                tar_start=1+lag2;
            else
                s1_start=1;
                s2_start=1;
                tar_start=1+lag2;
            end
            
            S1t = Data(s1_start:s1_start+nTuples-1,n1);
            S2t = Data(s2_start:s2_start+nTuples-1,n2); %same node, different lag
            Tar =  Data(tar_start:tar_start+nTuples-1,Target);
            
            
            if std(S1t)<lim && std(S2t)<lim && std(Tar)<lim
                continue
            end
            
            pdf = compute_pdfGUI([S1t S2t Tar],N,bin_scheme,...
                [Range(:,n1) Range(:,n2) Range(:,Target)],method, [z_effect(n1+1) z_effect(n2+1) z_effect(Target+1)]);
            info = compute_info_measures(pdf);
            
        R_allpairs(pairct,Target)=info.R;
        S_allpairs(pairct,Target)=info.S;
        I_ind1_allpairs(pairct,Target)=info.I_x1y;
        I_ind2_allpairs(pairct,Target)=info.I_x2y;
        TI_allpairs(pairct,Target) = info.T/info.Itot;
        Itot_allpairs(pairct,Target) = info.Itot;
            
        end 
        pairct=pairct+1;
    end
 
   timebar(timer,v1/nvars)   
    
end

%for each individual link, compute unique information (min val of I_dom - R)
for n =1:nvars    
    for Target=1:nvars
      search_ind1=[]; search_ind2 =[];
      search_ind1 = find(pairvars(:,1)==n);
      search_ind2 = find(pairvars(:,2)==n);
      
      buddy_1 = pairvars(search_ind1,2);
      buddy_2 = pairvars(search_ind2,1);
      buddyvect = [buddy_1' buddy_2'];
      
      
      I_vals1 = I_ind1_allpairs(search_ind1,Target);
      I_vals2 = I_ind2_allpairs(search_ind2,Target);
      R_vals1 = R_allpairs(search_ind1,Target);
      S_vals1 = S_allpairs(search_ind1,Target);
      R_vals2 = R_allpairs(search_ind2,Target);
      S_vals2 = S_allpairs(search_ind2,Target);
      Itot_vals1 = Itot_allpairs(search_ind1,Target);
      Itot_vals2 = Itot_allpairs(search_ind2,Target);
      TI_vals1 = TI_allpairs(search_ind1,Target);
      TI_vals2 = TI_allpairs(search_ind2,Target);
      
      Rvect = [R_vals1' R_vals2'];
      Svect = [S_vals1' S_vals2'];
      Itotvect = [Itot_vals1' Itot_vals2'];
      TIvect = [TI_vals1' TI_vals2'];
      TIvect(TIvect==0)=nan;
      
      [R(n,Target) R_index] = max(Rvect);
      [S(n,Target) S_index] = max(Svect);
      [TI(n,Target) TI_index] = min(TIvect);
      [Itot(n,Target) Itot_index] = max(Itotvect);
      
      S_pair(n,Target) = buddyvect(S_index);
      R_pair(n,Target) = buddyvect(R_index);
      TI_pair(n,Target) = buddyvect(TI_index);
      Itot_pair(n,Target) = buddyvect(Itot_index);
      
      U_vals1 = I_vals1-R_vals1;
      U_vals2 = I_vals2-R_vals2;
      
      U_vals1(U_vals1==0)=[]; U_vals2(U_vals2==0)=[];
      
      if length(U_vals1)>=1 || length(U_vals2)>=1
      U(n,Target) = min([U_vals1' U_vals2']);
      
      end
           
    end
       
end


S_pair(S==0)=0;
R_pair(R==0)=0;

TI(isnan(TI))=0;

close(timer)

%% save output products in entropy structure
entropy.H_x_1 =             H_x_1;              %entropy based on 1D pdf
entropy.H_x_2 =             H_x_2;              %entropy based on 2D pdf
entropy.I_tau =             I_lags;             %mutual info I(lag,source,target)
entropy.I_tau_normbyH =     I_normbyH;          %normalized mutual info I(lag,source,target)
entropy.I_dom =             I_dom;              %maximum lagged mutual information (bits)
entropy.I_dom_normbyH =     I_dom_normbyH;      %maximum lagged mutual information (normalized by H)
entropy.I_dom_lag=          I_tau;              %lag associated with I_dom (timesteps)
entropy.TE_tau=             TE;                 %transfer entropy TE(lag, source, target)
entropy.TE  =               TE_T;               %transfer entropy asociated with I_dom (TE at lag of strongest link)
entropy.TE_normbyItot =     TE_T_normbyItot;    %transfer entropy normalized by Itot = I(Source;Target_history;Target)

entropy.p_value =           I_pvalue;         %p-value of statistical significance of links

entropy.Itot =              Itot;             %maximum total info from 2 sources I(source,alt_source:target)
entropy.TI =                TI;               %minimum T/I = unique+synergy index TI(source,target)
entropy.U =                 U;                %unique information associated with dominant I links U(source,target)
entropy.S =                 S;                %max synergistic information associated with dominant I links S(source,target)
entropy.R =                 R;                %redundant information associated with dominant I links R(source,target) 

entropy.pairvars =          pairvars;  %source variable pairs

entropy.Itot_allpairs =     Itot_allpairs;    %maximum total information to each target from all source pairs
entropy.TI_allpairs =       TI_allpairs;      %T/I information to each target from all source pairs
entropy.I_ind1_allpairs =   I_ind1_allpairs;
entropy.I_ind2_allpairs =   I_ind2_allpairs;
entropy.S_allpairs =        S_allpairs;       %synergistic information to each target from all source pairs
entropy.R_allpairs =        R_allpairs;       %redundant information to each target from all source pairs 

entropy.S_pair =            S_pair;            %other source node providing maximum syn info SpairT(a,b)=c means a provides Syn info to b along with c
entropy.R_pair =            R_pair;            %other source node providing maximum red info RpairT(a,b)=c means a and c provide redundant info to b
entropy.Itot_pair =         Itot_pair;           %other source node associated with max total I Itot_pair(a,b)=c means a and c provide highest Ua+Uc+R+S to b
entropy.TI_pair =           TI_pair;           %other source node associated with minimum T/I TIpaitT(a,b)=c means a and c provide lowest (Ua+S)/Itot

entropy.N_sources =         N_sources;          %number of sources detected to each target (same node can transmit at mult lags)
entropy.N_sourcenodes =     N_sourcenodes;      %number of individual source nodes to each target
entropy.I_inst =            I_inst;             %lag zero mutual information
entropy.I_inst_normbyH=     I_inst_normbyH;     %normalized lag zero mutual information
entropy.I_inst_sigthresh =  I_inst_sig;   %statistical significance threshold for mutual information

end

