function [entropy] = EntropyFun_18GUI_zerolag(mi,Data,seg)
%Function to compute various entropy and transfer entropy measures
%uses KDE method
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
%September 30, 2015: updated compute_info_measures2 for output as structure
%January 18, 2016: updates for GUI - input h (smoothing parameter) to KDE
%computation, option for method in mi - fixed binning or KDE
%February 2016: update options: mi.netopt (1=full network, 2 = H and I
%only, 3 = H only (only entropy of nodes)
%February 2016: alter redundancy: 2 redundancy matrices: R_tau is
%redundancy between a single source to a target node at different time
%lags, R_T is redundancy between dominant link from each source to a target

%Update: altered version to consider zero-lag I as possible dominant link
%for T/I, other measures, depending on mi.ZeroLagOpt

nvars = mi.nvars;
nTests = mi.nTests;
N = mi.N;
lagvect = mi.lagvect;
bin_scheme=mi.bin_scheme;
Range = mi.Range;
maxS = mi.maxS;
method = mi.method;
opt = mi.netopt;
z_opt = mi.ZeroLagOpt;

if z_opt==1
lagvect =[0 lagvect];
end

h1D = mi.KDEparams.h1;
h2D = mi.KDEparams.h2;
h3D = mi.KDEparams.h3;

ndata = size(Data,1);
nlags = length(lagvect);
lim = eps*N^3; % machine limit of matlab * number of bins^3
N_sources = zeros(nvars,1);
N_sourcenodes = zeros(nvars,1);
TE = zeros(nlags,nvars,nvars);
TE_normbyItot = TE;
I_lags = zeros(nlags,nvars,nvars);
I_inst = zeros(nvars,nvars);
Hxtemp = zeros(nlags,nvars,nvars);
Hytemp = zeros(nlags,nvars,nvars);
H_x = zeros(nvars,1);

%initially set redundant and unique information to zero
R_tau = zeros(nlags,nvars,nvars);
S_tau = R_tau;
U_tau = R_tau;

S_T = zeros(nvars);
R_T = S_T;
U_T = S_T;
Spair_T = S_T;  %synergy buddies

S_T_normbyItot = S_T;
U_T_normbyItot = S_T;
R_T_normbyItot = S_T;
TE_T_normbyItot = S_T;

TI_T = zeros(nvars);
TE_T = TI_T;
Itot_T = zeros(nvars);
I_dom_strength = zeros(nvars);
I_dom_tau = zeros(nvars);
I_unique_tau = zeros(nvars);

%% flag variables with constant values or NaN values, don't use
flag = zeros(1,nvars);
for Trans=1:nvars
    X= Data(:,Trans);
    numNaNs = sum(isnan(X));
    num9999 = sum(X==-9999);
    %check for all zero or all constant values
    if std(X)<lim || numNaNs >0 || num9999 >0
        flag(Trans)=1;
    end
end

%% compute 1D entropy values for each node
for Trans = 1:nvars
    
    if flag(Trans)==1
        continue
    end
    X = Data(:,Trans);
    
    pdf = compute_pdfGUI(X,N,bin_scheme, Range(:,Trans),method,h1D(Trans));
    info = compute_info_measures2(pdf);
    H_x(Trans) = info.Hx;
end

if opt==1                                                   %done, exit function
    entropy.H_x = H_x;
    return
end

%% compute MI and lagged MI for all lags and pairs

timer=timebar(sprintf('MI for segment %d',seg),'Network Timer');
ct=0;
timebar(timer,ct/nvars^2);

for Trans = 1:nvars %transmitters
    
    X = Data(:,Trans);
    
    if flag(Trans)==1
        ct=ct+nvars;
        continue
    end
    
    for Rec = 1:nvars %receivers
        %fprintf('trans = %d, rec = %d\n',Trans,Rec)
        
        Y = Data(:,Rec);
        if flag(Rec)==1
            ct=ct+1;
            continue
        end
   
        for t=1:nlags %may or may not include zero lag           

            lagt = lagvect(t);

            nTuples = ndata-lagt-1;
            XtSTART=1;
            YfSTART=1+lagt;
            YwSTART = max(1,lagt);
            
            Xt=X(XtSTART:XtSTART+nTuples-1);        %Leading Node Xt (lag tau earlier than present)
            Yf=Y(YfSTART:YfSTART+nTuples-1);        %Led Node Yf (one timestep in future)
            Yw=Y(YwSTART:YwSTART+nTuples-1);        %current node Yw
            
            %first compute mutual information for all lags
            % fprintf('computing lagged I Trans %d Rec %d lag %d\n',Trans,Rec,lagt)
            if std(Xt)>lim && std(Yf)>lim
                pdf = compute_pdfGUI([Xt Yf],N,bin_scheme, [Range(:,Trans)  Range(:,Rec)],method,...
                    [h2D(Trans) h2D(Rec)]);
                info = compute_info_measures2(pdf);
                I_lags(t,Trans,Rec)=info.I;
                
                Hxtemp(t,Trans,Rec) = info.Hx;
                Hytemp(t,Trans,Rec) = info.Hy;   
                
            %shuffle to compute significance
            if nTests>0
                if t==1 %do shuffled significance testing for zero lag only
                    I_shuff=zeros(1,nTests);
                    for test=1:nTests
                        Yf_shuff = randsample(Yf,length(Yf));
                        pdfshuff = compute_pdfGUI([Xt Yf_shuff],N,bin_scheme, ...
                            [Range(:,Trans)  Range(:,Rec)],method,[h2D(Trans) h2D(Rec)]);
                        infoshuff = compute_info_measures2(pdfshuff);
                        I_shuff(test) = infoshuff.I;
                    end
                    I_shuff_sig = mean(I_shuff)+3*std(I_shuff);
                end
 
                I_lags(t,Trans,Rec) = max(info.I-I_shuff_sig,0);        %not sig = 0
            end
            
            % for significant lagged I, also compute Transfer Entropy TE
            if I_lags(t,Trans,Rec)>0  %&& t>1 %MAJID now computing zero lag TE
                pdf = compute_pdfGUI([Xt Yw Yf],N,bin_scheme,...
                    [Range(:,Trans) Range(:,Rec) Range(:,Rec)],method,...
                    [h3D(Trans) h3D(Rec) h3D(Rec)]);
                info = compute_info_measures2(pdf);
                TE(t,Trans,Rec)=info.T;
                if info.I_tot>0
                TE_normbyItot(t,Trans,Rec) = info.T/info.I_tot; %normalized by total information shared
                end
               
            end
            
            end
            
        end                                                         %lags
        
        if z_opt ==1 %set I_inst to first index of I_lags
            I_inst(Trans,Rec) = I_lags(1,Trans,Rec);
        else   %otherwise compute instantaneous mutual info
            pdf = compute_pdfGUI([X Y],N,bin_scheme, [Range(:,Trans)  Range(:,Rec)],method,...
                [h2D(Trans) h2D(Rec)]);
            info = compute_info_measures2(pdf);
            I_lags(t,Trans,Rec)=info.I;
            if nTests>0
                I_shuff=zeros(1,nTests);
                for test=1:nTests
                    Y_shuff = randsample(Y,length(Yf));
                    pdfshuff = compute_pdfGUI([X Y_shuff],N,bin_scheme, ...
                        [Range(:,Trans)  Range(:,Rec)],method,[h2D(Trans) h2D(Rec)]);
                    infoshuff = compute_info_measures2(pdfshuff);
                    I_shuff(test) = infoshuff.I;
                end
                I_shuff_sig = mean(I_shuff)+3*std(I_shuff);
                
                I_inst(Trans,Rec) = max(info.I-I_shuff_sig,0);        %not sig = 0
            end
            
        end
        
        ct=ct+1;
        timebar(timer,ct/nvars^2)
    end                                                             %receivers
end                                                                 %transmitters


%normalize lagged I by entropy of variables
I_normbyH = I_lags ./ min(Hxtemp,Hytemp);
I_normbyH(isnan(I_normbyH))=0;
I_normbyH(I_normbyH>1)=0;

close(timer)

if opt==2                                                   %done, exit function
    entropy.H_x = H_x;
    entropy.Hxtemp = Hxtemp;
    entropy.Hytemp = Hytemp;
    entropy.I_lags = I_lags;
    entropy.I_normbyH = I_normbyH;
    entropy.I_inst = I_inst;
    return
end

%% use lagged MI to compute synergies/redundancies between sources

timer=timebar(sprintf('Redundancy for segment %d',seg),'Network Timer');
timebar(timer,0/nvars);

for Rec = 1:nvars                                       %cycle through all receiving (target) nodes
    
    if flag(Rec)==1 || sum(sum(I_lags(:,:,Rec)))==0
        continue                                        %go to next target node if Rec has no sources
    end

    %Determine top (maxS) source nodes for each target, and time lags
    matrix = reshape(I_lags(:,:,Rec),nlags,nvars);
    
    sig_Sources=[]; sig_lags=[]; sig_Ivalues=[];
    [sig_lags, sig_Sources]=find(matrix>0);             %sig_lags are indices of lags, not actual lags
    sig_Ivalues = matrix(matrix>0);
    
    [sig_I_sort, ind]=sort(sig_Ivalues,'descend');
    if length(sig_Ivalues)> maxS                        %cut off smaller links depending on maxS
        ind = ind(1:maxS);
        sig_I_sort = sig_I_sort(1:maxS);
        sig_lags=sig_lags(ind);                         %shortened vector of source lags
        sig_Sources=sig_Sources(ind);                   %shortened vector of source nodes
    end
    
    numsources = length(sig_Sources);                   %number of sources nodes <= maxS
    %fprintf('TARGET node %d has %d sources\n',Rec,numsources)
    N_sources(Rec)=numsources;
    N_sourcenodes(Rec)=length(unique(sig_Sources));
    
    if numsources==1                                    %target only has 1 source node/lag pair: all unique
        U_tau(sig_lags,sig_Sources,Rec)=I_lags(sig_lags,sig_Sources,Rec);
        TI(sig_lags,sig_Sources,Rec)=1;
        continue
    end
    
    source_nodes = unique(sig_Sources);                     %number of distinct source nodes
    redflag = zeros(nvars,length(lagvect));               %redundancy flag (to skip redundant nodes)
    
    %Source Redundancy (redudancy for each source node)
    for j = 1:length(source_nodes)                          %for each source node, consider each lag time
        
        node = source_nodes(j);
        lags = sig_lags(sig_Sources==node);
        if length(lags) == 1
            continue                                        %only one lag for that source, move to next source node
        end
        
        
        pairs = combnk(lags,2);  %nx2 matrix of all possible (lag1,lag2) pairs
        
        for p = 1:size(pairs,1)
            
            lag1 = max(pairs(p,:));
            lag2 = min(pairs(p,:));
            lag1t = lagvect(lag1); %actual lag times in lagvect
            lag2t = lagvect(lag2);
            
            if redflag(node,lag1)==1 || redflag(node,lag2)==1
                continue %move to next pair
            end

            nTuples = ndata-lag1t-1;
            
            S1START=1;
            S2START=1+lag1t-lag2t;
            YfSTART=1+lag1t;
            
            S1t = Data(S1START:S1START+nTuples-1,node);
            S2t = Data(S2START:S2START+nTuples-1,node); %same node, different lag
            Yf =  Data(YfSTART:YfSTART+nTuples-1,Rec);
            
            if std(S1t)<lim && std(S2t)<lim && std(Yf)<lim
                continue
            end
            
            pdf = compute_pdfGUI([S1t S2t Yf],N,bin_scheme,...
                [Range(:,node) Range(:,node) Range(:,Rec)],method,...
                [h3D(node) h3D(node) h3D(Rec)]); 
            info = compute_info_measures2(pdf);
           
            %Redundancy
            R_tau(lag1,node,Rec)=info.R;
            R_tau(lag2,node,Rec)=info.R;
                  
            if info.U1 <=0.01 * info.R            %if unique information is very small compared to redundancy
                redflag(node,lag1)=1;             %classify this (source, lag) pair as redundant
                R_tau(lag1,node,Rec)=info.I_lags;  %assume all shared info is redundant
                I_lags(lag1,node,Rec)=info.I_lags;  %re-define I_lags
                U_tau(lag1,node,Rec)=info.U1;
            elseif info.U1>U_tau(lag1,node,Rec)
                U_tau(lag1,node,Rec)=info.U1;
            end
            
            if info.U2 <=0.01 * info.R            %if unique information is very small compared to redundancy
                redflag(node,lag1)=1;             %classify this (source, lag) pair as redundant
                R_tau(lag2,node,Rec)=info.I_x2y;  %assume all shared info is redundant
                I_lags(lag2,node,Rec)=info.I_x2y;  %re-define I_lags
                U_tau(lag1,node,Rec)=info.U1;
            elseif info.U2>U_tau(lag2,node,Rec)
                U_tau(lag2,node,Rec)=info.U2;
            end
            
            if info.S > S_tau(lag1,node,Rec)      %maximum S for that lag
                S_tau(lag1,node,Rec) = info.S;
            end
            
            if info.S > S_tau(lag2,node,Rec)
                S_tau(lag2,node,Rec) = info.S;
            end
            
        end                                                 %end source pairs 
    end                                                     %end unique source nodes   
end                                                         %end Receivers

%matrix of dominant links (a) in strength and (b) in uniqueness
for Rec =1:nvars
    for Trans = 1:nvars
        vect = reshape(I_normbyH(2:end,Trans,Rec),1,nlags);
        sprintf('%4.3f\n',max(vect));
        if max(vect)>0
            [I_dom_strength(Trans,Rec), ind_lag] = max(vect);    %maximum strength
            I_dom_tau(Trans,Rec) = lagvect(ind_lag(1));
            TE_T(Trans,Rec) = TE(ind_lag(1)+1,Trans,Rec);
            TE_T_normbyItot(Trans,Rec) = TE_normbyItot(ind_lag(1)+1,Trans,Rec);
            vect2 = reshape(U_tau(:,Trans,Rec),1,nlags);
            [~, ind_lag] = max(vect2); %maximum uniqueness
            I_unique_tau(Trans,Rec)=lagvect(ind_lag(1));
        end
    end
end

I_dom = I_dom_strength;
Lag_dom = I_dom_tau;
  
for Rec = 1:nvars                                           %cycle back through recievers for network R,S,U
       
    %For each dominant source (lag, node) pair, compute U, R, S Info
    US1=[]; US2 =[]; SynS1S2=[]; RedS1S2=[];I_lagged_S1=[]; I_lagged_S2=[];
    TN=[];
    
    ct=1;
    %pairwise computation of redundancy, synergy
    for Trans1 = 1:(nvars-1)
        for Trans2 = (Trans1+1):nvars
            
            if I_dom(Trans1,Rec)==0 || I_dom(Trans2,Rec)==0
                continue
            end
            
            lagS1 = Lag_dom(Trans1,Rec);
            lagS2 = Lag_dom(Trans2,Rec);

            nTuples = ndata-max(lagS1,lagS2)-1;
            
            S2START=1; S1START=1;
            if lagS1>lagS2
                S2START=1+lagS1-lagS2;
            elseif lagS1<lagS2
                S1START=1+lagS2-lagS1;
            end
            YfSTART = 1+max(lagS1,lagS2);
            
            %% Information Measures
            S1t=Data(S1START:S1START+nTuples-1,Trans1);
            S2t=Data(S2START:S2START+nTuples-1,Trans2);
            Yf=Data(YfSTART:YfSTART+nTuples-1,Rec);
            
            pdf = compute_pdfGUI([S1t S2t Yf],N,bin_scheme,...
                [Range(:,Trans1) Range(:,Trans2) Range(:,Rec)],method,...
                [h3D(Trans1) h3D(Trans2) h3D(Rec)]);
            
            info = compute_info_measures2(pdf);
            I_lagged_S1(ct) = info.I_lags;
            I_lagged_S2(ct) = info.I_x2y;
            I_totS1S2(ct) = info.I_tot;
            RedS1S2(ct) = info.R;
            SynS1S2(ct) = info.S;
            US1(ct) = info.U1;
            US2(ct) = info.U2;
            dI(ct)=info.dI;
            T(ct) = info.T;
            
            ToverTOTI(ct) = T(ct)/I_totS1S2(ct);
            
            TN(ct,:)=[Trans1 lagS1 Trans2 lagS2];
            ct=ct+1;
        end
    end
     
    %for each source, Unique, Redundant, and Synergistic Info
    for Trans =1:nvars

        if I_dom(Trans,Rec)==0 || size(TN,1)==0
            continue
        end
        lag = Lag_dom(Trans,Rec);
        
        indA=[]; indB=[];
        indA = find(TN(:,1)==Trans & TN(:,2)==lag);
        indB = find(TN(:,3)==Trans & TN(:,4)==lag);
        
        TOTIAB =[I_totS1S2(indA) I_totS1S2(indB)];
        TAB = [ToverTOTI(indA) ToverTOTI(indB)];
        UAB = [US1(indA) US2(indB)];
        SAB =[SynS1S2(indA) SynS1S2(indB)];

        SynBuddy = [TN(indA,3); TN(indB,1)]; %other node providing synergy with Trans
        RAB =[RedS1S2(indA) RedS1S2(indB)];
        
        [R_T(Trans,Rec), ind] = max(RAB); %index for max redundancy
        if TOTIAB(ind)>0
            R_T_normbyItot(Trans,Rec)=max(RAB)./TOTIAB(ind);
            
            
            if UAB(ind)>.01*max(RAB) %if unique information is significant compared to redundancy
                U_T(Trans,Rec)= UAB(ind);
                U_T_normbyItot(Trans,Rec) = UAB(ind)./TOTIAB(ind);
            end
        end
        
        [S_T(Trans,Rec), ind]= max(SAB); %new index for maximum synergy
        Spair_T(Trans,Rec)=SynBuddy(ind); %node that accompanies Trans to provide Syn Info
        
        if TOTIAB(ind)>0
        S_T_normbyItot(Trans,Rec)=max(SAB)./TOTIAB(ind);
        end
        
        TI_T(Trans,Rec)=min(TAB);
        Itot_T(Trans,Rec)=max(TOTIAB);
    end
    
    timebar(timer,Rec/nvars)
end                                                                         %end of cycle

close(timer)

%% save output products in entropy structure
entropy.H_x =           H_x;
entropy.I_lags =         I_lags;
entropy.I_normbyH =     I_normbyH;
entropy.Itot_T =        Itot_T;
entropy.TI_T =          TI_T;
entropy.U_tau =         U_tau;
entropy.U_T =           U_T;
entropy.U_T_normbyItot = U_T_normbyItot;
entropy.S_tau =         S_tau;
entropy.S_T =           S_T;
entropy.Spair_T =       Spair_T;
entropy.S_T_normbyItot = S_T_normbyItot;
entropy.R_tau =         R_tau;
entropy.R_T =           R_T;
entropy.R_T_normbyItot = R_T_normbyItot;
entropy.I_dom =         I_dom_strength;
entropy.I_dom_tau=      I_dom_tau;
entropy.I_unique_tau =  I_unique_tau;
entropy.TE=             TE;
entropy.TE_T =          TE_T;
entropy.TE_T_normbyItot = TE_T_normbyItot;
entropy.N_sources =     N_sources;
entropy.N_sourcenodes = N_sourcenodes;
entropy.I_inst =        I_inst;

end

