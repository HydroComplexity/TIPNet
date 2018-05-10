function AllStats = TotalStatsFunctionGUI(netdata,mi,entres)
%TotalStatsFunction
%compute network statistics for a network

%netdata: network data structure (one cell for each network)
%mi: model information structure
%entres: entropy results structure (one cell for each network)

%Updated 9/10/15 by Allison Goodwell
%Updated 2/4/16 for GUI
%Updated 2/17/16 for updated EntropyFun_18GUI function

lagvect= mi.lagvect;
n_nets_total = size(netdata,2);
nvars = mi.nvars;
nlags = length(lagvect);

I_mut = zeros(nvars,nvars,n_nets_total);

%% Statistics over all networks

if mi.netopt==2
  
    
    
elseif mi.netopt==3
    

for n = 1:n_nets_total
 
        I_max(:,:,n)=entres{n}.I_dom;
        I_lag(:,:,n)=entres{n}.I_dom_lag;
        U_Imax(:,:,n)=entres{n}.U;
        R_Imax(:,:,n)=entres{n}.R;
        TI_Imax(:,:,n)=entres{n}.TI;
        TE_Imax(:,:,n)=entres{n}.TE;
        S_Imax(:,:,n) = entres{n}.S;       
        I_mut(:,:,n)=entres{n}.I_inst;
        I_mut_normbyH(:,:,n) = entres{n}.I_inst_normbyH;
        
                          
end

for i =1:nvars
    for j = 1:nvars
        Imax_var = reshape(I_max(i,j,:),1,n_nets_total);
        ind=[];
        ind = find(Imax_var>0 & ~isnan(Imax_var));
        AllStats.Link_freq(i,j) = n_nets_total - length(ind);
        
        %all-network averages
        AllStats.I_inst(i,j) = mean(I_mut(i,j,:));                  %mean inst mutual information
        AllStats.I_inst_normbyH(i,j) = mean(I_mut_normbyH(i,j,:));  %mean normalized inst mutual info
        AllStats.I_normbyH(i,j) = mean(I_max(i,j,ind));             %mean dominant lagged I
        AllStats.I_lag(i,j) = mean(I_lag(i,j,ind));                 %mean lag associated with dominant I
        AllStats.U(i,j) = mean(U_Imax(i,j,ind));                    %mean U associated with dominant I
        AllStats.TI(i,j) = mean(TI_Imax(i,j,ind));                  %mean TI associated with dominant I
        AllStats.TE(i,j) = mean(TE_Imax(i,j,ind));                  %mean TE associated with dominant I
        AllStats.R(i,j) = mean(R_Imax(i,j,ind));                    %mean R associated with dominant I
        AllStats.S(i,j) = mean(S_Imax(i,j,ind));                    %mean S associated with dominant I          
    end
end


I_sum = zeros(nlags,nvars,nvars);
TE_sum = zeros(nlags,nvars,nvars);

I_normbyH_t = zeros(nlags,nvars,nvars);
TE_t  = zeros(nlags,nvars,nvars);

for t=1:nlags
    
    for n = 1:n_nets_total
        
        I_sum(t,:,:)=I_sum(t,:,:)+entres{n}.I_tau_normbyH(t,:,:);
        TE_sum(t,:,:)=TE_sum(t,:,:)+entres{n}.TE_tau(t,:,:);  
    end
    
    I_normbyH_t(t,:,:)=I_sum(t,:,:)./n_nets_total;
    TE_t(t,:,:)=TE_sum(t,:,:)./n_nets_total;
   
end

AllStats.I_normbyH_t=I_normbyH_t;               %for each lag: mean lagged I
AllStats.TE_t=TE_t;                             %for each lag: mean TE

end

