function [data_out] = ProcessFunction(data,rm_outliers,norm)
%ProcessFunction: to pre-process data for network measures
%data = dataset
%rm_outliers: 1 to remove outliers
%norm: 1 to normalize between 0 and 1

if rm_outliers ==1
    
    upbound = 75; %higher to only get rid of most extreme outliers
    lowbound = 25;
    
    for j = 1:size(data,2)
        
        
        dat = data(:,j);
        dat2 = round(dat.*100)./100;
        m = mode(dat2);
        if length(dat2(dat2==m))>.1*length(dat)
        dat_nozeros = dat;
        dat_nozeros(dat2==m)=[];
        else
        dat_nozeros = dat;
        end
        
        p75 = prctile(dat_nozeros,upbound);
        p25 = prctile(dat_nozeros,lowbound);
        IQR = p75-p25;
        
        uplim = p75+1.5.*IQR;
        lowlim = p25-1.5.*IQR;
        
        dat(dat>uplim)=uplim;
        dat(dat<lowlim)=lowlim;
        
        
        
        data_out(:,j)=dat;
    end
else
    
    data_out = data;
end

if norm ==1
    
    minval = min(data_out);
    maxval = max(data_out);
    for j = 1:size(data_out,2)
        if maxval(j)-minval(j)>0 %if range >0
        data_out(:,j) = (data_out(:,j) - minval(j))./(maxval(j)-minval(j));
        end
    end
    
    
end



