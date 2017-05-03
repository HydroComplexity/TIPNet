function [Data] = chaotic_generate(ndata,lagmat,epsilon,epsilonz)

% generate X and Y data with logistic map
% ndata = number of points to generate
% adjmat = adjacency matrix with (1 for link, 0 for no link)
% lagmat = matrix with coupling lags
% epsilon = coupling strength (uniform for all links) between 0 and 1
% epsmat = only used in synchronized systems case (11) - matrix of coupling
% strengths (bi-modal)

% modified 9/18/14 to add coupling strength (epsilon)

trash = 1000; 
a = 4;
nvars = size(lagmat,1);
Data = rand(ndata+trash,nvars); 
start = max(2,max(max(lagmat))+1);

for t = start:trash+ndata+(start-1)
    for i = 1:nvars %i is receiving from j
        sum_strong =0;
        k_strong =0; %in-degree of node
        for j = 1:nvars %j is transmitting to i        
            if lagmat(i,j) > 0  %if connected
                sum_strong = sum_strong + a*Data(t-lagmat(i,j),j)*(1-Data(t-lagmat(i,j),j));
                k_strong =k_strong+1; 
            end          
        end
      
        term_ind = (1-epsilon)*a*Data(t-1,i)*(1-Data(t-1,i)); %independent term
        term_link = epsilon*(1-epsilonz)*sum_strong/max(k_strong,1); %linked term
        term_noise = epsilon*epsilonz*rand();
        
        Data(t,i)=term_ind + term_link + term_noise;
      
    end
end
    
    %clear first values
    if trash>0
    Data(1:start-1+trash,:)=[];
    end
    
 end
% end of function
