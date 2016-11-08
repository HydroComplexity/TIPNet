function [pdf, Coords]= compute_pdfGUI(Data,N,bin_scheme, Range,method,h)
% Compute pdf 1D, 2D, or 3D
% use KDE method, Epanichov (Silverman 1986) kernel to compute 1d, 2d, 3d
% pdfs

%N = number of grid points at which to compute kernel
%bin_scheme: local or global
%Range: used for global binning scheme
%Updated 2/19/15 to normalize all measures (Hnorm added)
%Updated 6/15/15 to split functions between finding pdf and calculating
%information measures
%1/18/16: switch h (smoothing param) accoring to GUI input: h is 1xDim
%matrix
%1/18/16: working on correction for zero effect - done
%1/25/16: add choice of method for pdf: KDE or fixed binning (i.e. for
%binary)
pdf=0;

nTup = size(Data,1);
dim = size(Data,2);

%measures to determine at which grid points to consider contributions from
%each data point: delta, h, xo

Coords = zeros(dim,N);
Edges = zeros(dim,N+1);

%compute non-zero fraction for each dimension to determine whether to
%correct for zero-effect (if more than 10% zero values)
atom = zeros(1,dim);
for i =1:dim
    kx(i) = sum(Data(:,i)~=0)./nTup;
    h(i)=h(i)./sqrt(kx(i));
end

atom(kx<.9)=1;
kx(kx>=0.9)=1;

if sum(atom)>=1 && N < 10
    method = 'fixed'; %switched to fixed bins if many zero values or few bins
end

if sum(kx)==0 %%%%% all zero values sent (don't compute anything)
    
    if dim==1
        pdf = zeros(1,N);
        pdf(1)=1;
    elseif dim==2
        pdf=zeros(N,N);
        pdf(1,1)=1;
    elseif dim==3
        pdf=zeros(N,N,N);
        pdf(1,1,1)=1;
    end
    
    return %%%%%% exit function
end


for i = 1:dim
        if atom(i)==0 %no atom at zero, proceed as usual
            if strcmp(bin_scheme,'local')==1
                xo = min(Data);                                 %xo(dim)
                Edges(i,:)=linspace(min(Data(:,i)), max(Data(:,i)),N+1);
            elseif strcmp(bin_scheme,'global')==1
                xo = Range(1,:) ;                                 %xo(dim)
                Edges(i,:) = linspace(Range(1,i),Range(2,i),N+1);
            end
            Coords(i,:)=(Edges(i,1:end-1)+Edges(i,2:end))./2;
            
        else %atom at zero: first bin is only for zero values, evenly space other bins
            xo = zeros(size(min(Data)));
            if strcmp(bin_scheme,'local')==1
                Edges(i,:)=[0 linspace(min(Data(:,i)+10^-8), max(Data(:,i)),N)];
            elseif strcmp(bin_scheme,'global')==1
                Edges(i,:) = [0 linspace(Range(1,i)+10^-8,Range(2,i),N)];
            end
            Coords(i,:)=(Edges(i,1:end-1)+Edges(i,2:end))./2;
            Coords(i,1)=0;
            

        end
end


delta = (Coords(:,end)-Coords(:,end-1))' ;        %delta(dim)

if strcmp(method,'KDE')==1 %%%%%%%%%%% KDE Method  %%%%%%%%%%%%%%%%%
    
    pdfcenter = zeros(1,N^dim);
    pdfx = zeros(1,N);
    pdfy = zeros(1,N);
    pdfz = zeros(1,N);
    
    pdfxy = zeros(1,N^2); %for 3D pdfs only
    pdfxz = zeros(1,N^2);
    pdfyz = zeros(1,N^2);
    
    for n =1:nTup %look at each data point separately
        
        dat = Data(n,:);
        
        %determine i,j,k indices at which to compute pdf values
            minind = floor((dat - xo - h)./ delta) -1;
            maxind = ceil((dat -xo + h)./ delta) +1;
            minind(minind<1 | isnan(minind))=1;
            maxind(maxind>N | isnan(maxind))=N;

        %if any index has zero effect, do not look at 1st coord
        if sum(atom)>=1
            minind=max(2,minind);
        end

        
        ICoords = Coords(1,minind(1):maxind(1));
        if dim >1
            JCoords = Coords(2,minind(2):maxind(2));
            if dim == 3
                KCoords = Coords(3,minind(3):maxind(3));
            end
        end
        
        %use c mex functions to calculate pdf
        if dim==1
            pdfadd = mdKDE_1d(dat,minind,maxind,ICoords,N,h);

            
            pdfadd(isnan(pdfadd))=0;
            pdfcenter = pdfcenter+pdfadd;
            
        elseif dim ==2
            
            if sum(atom)>0 && sum(dat~=0)==2 %if atom exists, p(x,y) for nonzero pts only
                pdfadd = mdKDE_2d(dat,minind,maxind,ICoords,JCoords,N,h);
            elseif sum(atom)==0 %if no atom, evaluate for all points
                pdfadd = mdKDE_2d(dat,minind,maxind,ICoords,JCoords,N,h);
            else
                pdfadd=0;
            end
            
            pdfadd(isnan(pdfadd))=0;
            pdfcenter = pdfcenter+pdfadd;
            %marginal pdfs: pxo(x) for x~=0,y=0 and pyo(y) for x=0,y~=0
            
            if atom(1)==1 && dat(1)==0 && dat(2)~=0 %pyo: zero effect for x
                pdfyadd = mdKDE_1d(dat(2),minind(2),maxind(2),JCoords,N,h(2));
                pdfyadd(isnan(pdfyadd))=0;
                pdfy = pdfy + pdfyadd;
            end
            
            if atom(2)==1 && dat(1)~=0 && dat(2)==0 %pxo: zero effect for y
                pdfxadd = mdKDE_1d(dat(1),minind(1),maxind(1),ICoords,N,h(1));
                pdfxadd(isnan(pdfxadd))=0;
                pdfx = pdfx + pdfxadd;
            end
            
        elseif dim ==3
            
            if sum(atom)>0 && sum(dat~=0)==3 %if atom, p(x,y,z) for nonzero pts only
                pdfadd= mdKDE_3d(dat,minind,maxind,ICoords,JCoords,KCoords,N,h);
                pdfadd(isnan(pdfadd))=0;
                pdfcenter = pdfcenter+pdfadd;
            elseif sum(atom)==0 %if no atom, p(x,y,z) for any points
                pdfadd = mdKDE_3d(dat,minind,maxind,ICoords,JCoords,KCoords,N,h);
                pdfadd(isnan(pdfadd))=0;
                pdfcenter = pdfcenter+pdfadd;
            elseif atom(1)==1 && dat(1)==0 && dat(2)~=0 && dat(3)~=0 %p(y,z)
                pdfyzadd = mdKDE_2d([dat(2) dat(3)],[minind(2) minind(3)],...
                    [maxind(2) maxind(3)],JCoords,KCoords,N,[h(2) h(3)]);
                pdfyzadd(isnan(pdfyzadd))=0;
                pdfyz = pdfyz+pdfyzadd;
            elseif atom(2)==1 && dat(1)~=0 && dat(2)==0  && dat(3)~=0 %p(x,z)
                pdfxzadd = mdKDE_2d([dat(1) dat(3)],[minind(1) minind(3)],...
                    [maxind(1) maxind(3)],ICoords,KCoords,N,[h(1) h(3)]);
                pdfxzadd(isnan(pdfxzadd))=0;
                pdfxz = pdfxz+pdfxzadd;
            elseif atom(3)==1 && dat(1)~=0 && dat(2)~=0  && dat(3)==0 %p(x,y)
                pdfxyadd = mdKDE_2d([dat(1) dat(2)],[minind(1) minind(2)],...
                    [maxind(1) maxind(2)],ICoords,JCoords,N,[h(1) h(2)]);
                pdfxyadd(isnan(pdfxyadd))=0;
                pdfxy = pdfxy+pdfxyadd;
            elseif atom(1)==1 && atom(2)==1 && dat(1)==0 && dat(2)==0 && dat(3)~=0 %p(z)
                pdfzadd = mdKDE_1d(dat(3),minind(3),maxind(3),KCoords,N,h(3));
                pdfzadd(isnan(pdfzadd))=0;
                pdfz = pdfz + pdfzadd;
            elseif atom(1)==1 && atom(3)==1 && dat(1)==0 && dat(2)~=0 && dat(3)==0 %p(y)
                pdfyadd = mdKDE_1d(dat(2),minind(2),maxind(2),JCoords,N,h(2));
                pdfyadd(isnan(pdfyadd))=0;
                pdfy = pdfy + pdfyadd;
            elseif  atom(2)==1 && atom(3)==1 && dat(1)~=0 && dat(2)==0 && dat(3)==0 %p(x)
                pdfxadd = mdKDE_1d(dat(1),minind(1),maxind(1),ICoords,N,h(1));
                pdfxadd(isnan(pdfxadd))=0;
                pdfx = pdfx + pdfxadd;
            end
            
        end
        
    end
    
    pdfcenter = pdfcenter./sum(sum(sum(pdfcenter))); %sum to 1
    pdfcenter = pdfcenter.*prod(kx); %area of center pdf = (kx*ky*kz)
    
    pdfcenter(isnan(pdfcenter))=0;
    
    if dim==2
        pdfy = (pdfy./sum(pdfy)).*((1-kx(1))*kx(2));
        pdfy(isnan(pdfy))=0;
        pdfx = (pdfx./sum(pdfx)).*((1-kx(2))*kx(1));
        pdfx(isnan(pdfx))=0;
        pdfcenter = reshape(pdfcenter,N,N);
        
    elseif dim==3
        pdfx = (pdfx./sum(pdfx)).*(kx(1)*(1-kx(2))*(1-kx(3)));
        pdfy = (pdfy./sum(pdfy)).*((1-kx(1))*kx(2)*(1-kx(3)));
        pdfz = (pdfz./sum(pdfz)).*((1-kx(1))*(1-kx(2))*kx(3));
        pdfxy = (pdfxy./sum(pdfxy)).*(kx(1)*kx(2)*(1-kx(3)));
        pdfxz = (pdfxz./sum(pdfxz)).*(kx(1)*(1-kx(2))*kx(3));
        pdfyz = (pdfyz./sum(pdfyz)).*((1-kx(1))*kx(2)*kx(3));
        
        pdfx(isnan(pdfx))=0;
        pdfy(isnan(pdfy))=0;
        pdfz(isnan(pdfz))=0;
        pdfxy(isnan(pdfxy))=0;
        pdfxz(isnan(pdfxz))=0;
        pdfyz(isnan(pdfyz))=0;
        
        
        pdfxy = reshape(pdfxy,N,N);
        pdfxz = reshape(pdfxz,N,N);
        pdfyz = reshape(pdfyz,N,N);
        
        pdfcenter = reshape(pdfcenter,N,N,N);
    end
    
    
    pdf = pdfcenter;
    
    if dim==1 && sum(atom)>0 %merge the marginal and central pdfs together
        pdf(1)=1-kx;
    elseif dim==2 && sum(atom)>0 %pdf(x,y)
        
        if atom(1)==1
            pdf(1,:)=pdfy;
        end
        if atom(2)==1
            pdf(:,1)=pdfx;
        end
        
        pdf(1,1)=(1-kx(1))*(1-kx(2));
        
    elseif dim==3 && sum(atom)>0
        
        if sum(atom)>0
            pdf(1,1,1)=(1-kx(1))*(1-kx(2))*(1-kx(3)); %p(x=0,y=0,z=0)
        end
        
        if atom(1)==1  %p(y,z)
            pdf(1,2:end,2:end)=pdfyz(2:end,2:end);
            
            if atom(2)==1
                pdf(1,1,2:end)=pdfz(2:end); %p(z)
            end
            if atom(3)==1
                pdf(1,2:end,1)= pdfy(2:end); %p(y)
            end
            
        end
        
        if atom(2)==1   %p(x,z)
            pdf(2:end,1,2:end)=pdfxz(2:end,2:end);
            if atom(3)==1
                pdf(2:end,1,1)= pdfx(2:end); %p(x)
            end
        end
        
        if atom(3)==1 %p(x,y)
            pdf(2:end,2:end,1)=pdfxy(2:end,2:end);
        end
        
    end
    
elseif strcmp(method,'fixed')==1 %%%%%%%%%%% Fixed Bin Method  %%%%%%%%%%%%%%%%%
    if dim==1
        pdfcenter=zeros(1,N);
    elseif dim==2
        pdfcenter=zeros(N,N);
    elseif dim==3
        pdfcenter=zeros(N,N,N);
    end
    
    pdfx = zeros(1,N);
    pdfy = zeros(1,N);
    pdfz = zeros(1,N);
    pdfxy = zeros(1,N^2); %for 3D pdfs only
    pdfxz = zeros(1,N^2);
    pdfyz = zeros(1,N^2);
 
    for i = 1:dim
        
        dat = Data(:,i);
        bindata = ones(size(dat));
        edges = Edges(i,:);
        
        for e = 1:N
            bindata(dat>=edges(e) & dat<edges(e+1))= e;
            if e==N
            bindata(dat>=edges(e+1))=e;  
            end
        end

        BinData(:,i)=bindata; %(i,j,k) bin numbers for each data point
    end
    
    
    
    if dim==1
        C=zeros(1,N);
        for n = 1:nTup
            dat = BinData(n);
            C(dat)=C(dat)+1;
        end
    elseif dim==2
        C=zeros(N,N);
        for n = 1:nTup
            dat = BinData(n,:);
            C(dat(1),dat(2))=C(dat(1),dat(2))+1;
        end
    elseif dim==3
        C=zeros(N,N,N);
        for n = 1:nTup
            dat = BinData(n,:);
            C(dat(1),dat(2),dat(3))=C(dat(1),dat(2),dat(3))+1;
        end
    end
 
pdf = C./nTup;  
    
end


end

