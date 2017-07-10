function [pdf, Coords]= compute_pdfGUI(Data,N,bin_scheme, Range,method)
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
%3/24/17: alter to use direct count for zero-effect instead of ind. assumption
%3/24/17: also alter to have "constant effect" instead of "zero effect"

nTup = size(Data,1);
dim = size(Data,2);

if dim==1
    pdf = zeros(1,N);
elseif dim==2
    pdf=zeros(N,N);
elseif dim==3
    pdf=zeros(N,N,N);
end

%%%%%%%%%%%%%%%%%%% Determine KDE smoothing parameters %%%%%%%%%%%%%%%%
for i =1:dim  
    dat = Data(:,i);
    Ra = range(dat);
    varterm = var(dat);
    varterm(isnan(varterm))=0;    
h1D(i)= 1.06 * nTup^(-1/5) .* varterm./Ra; 
h2D(i)= 1.77 * nTup^(-1/6) .* varterm./Ra; 
h3D(i) = 2.78 * nTup^(-1/7) .* varterm./Ra; 
end

%measures to determine at which grid points to consider contributions from
%each data point: delta, h, xo

Coords = zeros(dim,N);
Edges = zeros(dim,N+1);

%%%%%%%%%%%%%%%%%%%%% Determine bin coordinates %%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:dim

            xo = zeros(size(min(Data)));
            if strcmp(bin_scheme,'local')==1
                Edges(i,:)=[0 linspace(min(Data(:,i)+10^-8), max(Data(:,i)),N)];
            elseif strcmp(bin_scheme,'global')==1
                Edges(i,:) = [0 linspace(Range(1,i)+10^-8,Range(2,i),N)];
            end
            Coords(i,:)=(Edges(i,1:end-1)+Edges(i,2:end))./2;
            Coords(i,1)=0;
            

            

end

delta = (Coords(:,end)-Coords(:,end-1))' ;        %delta(dim)

%%%%%%%%%%%%%%%%% zero-effect accounting %%%%%%%%%%%%%%%%%%%%%%%%%

 for i =1:dim
     zval(i) = mode(Data(:,i)); %zero or constant value
     vect = Data(:,i);
     vectdiff = abs(vect-zval(i))./delta(i);
     vect(vectdiff<2)=zval(i);
     Data(:,i)=vect;
     
     kx(i) = sum(Data(:,i)~=zval(i))./nTup; %proportions of non-zero values
 end
 


if dim ==1
    fx = sum(Data~=zval(1))./nTup;
    f0 = sum(Data==zval(1))./nTup;
    h=h1D;
elseif dim==2
    f0 = sum(Data(:,1)==zval(1) & Data(:,2)==zval(2))./nTup;
    fx = sum(Data(:,1)~=zval(1) & Data(:,2)==zval(2))./nTup;
    fy = sum(Data(:,1)==zval(1) & Data(:,2)~=zval(2))./nTup;
    fxy = sum(Data(:,1)~=zval(1) & Data(:,2)~=zval(2))./nTup;
    h=h2D;
elseif dim==3
    f0 = sum(Data(:,1)==zval(1) & Data(:,2)==zval(2) & Data(:,3)==zval(3))./nTup;
    fx = sum(Data(:,1)~=zval(1) & Data(:,2)==zval(2) & Data(:,3)==zval(3))./nTup;
    fy = sum(Data(:,1)==zval(1) & Data(:,2)~=zval(2) & Data(:,3)==zval(3))./nTup;
    fz = sum(Data(:,1)==zval(1) & Data(:,2)==zval(2) & Data(:,3)~=zval(3))./nTup;
    fxy = sum(Data(:,1)~=zval(1) & Data(:,2)~=zval(2) & Data(:,3)==zval(3))./nTup;
    fxz = sum(Data(:,1)~=zval(1) & Data(:,2)==zval(2) & Data(:,3)~=zval(3))./nTup;
    fyz = sum(Data(:,1)==zval(1) & Data(:,2)~=zval(2) & Data(:,3)~=zval(3))./nTup;
    fxyz = sum(Data(:,1)~=zval(1) & Data(:,2)~=zval(2) & Data(:,3)~=zval(3))./nTup;
    h=h3D;
end

for i =1:dim
            diff =abs(Coords(i,:)-zval(i));
            mergecoord(i) = find(diff==min(diff));
end


if sum(kx)==0 %%%%% all zero or const values sent (don't compute anything)    
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



%%%%%%%%%%%%%%%%%%%%%%%% KDE Method  %%%%%%%%%%%%%%%%%%%%%%%
if strcmp(method,'KDE')==1 
    
    pdfcenter = zeros(1,N^dim); %3D pdf (no zero values)
    pdfx = zeros(1,N);          %1D pdf (line where only x>0)
    pdfy = zeros(1,N);          %1D pdf (line where only y>0)          
    pdfz = zeros(1,N);          %1D pdf (line where only z>0)
    
    pdfxy = zeros(1,N^2);       %2D pdf (wall where x>0,y>0)
    pdfxz = zeros(1,N^2);       %2D pdf (wall where x>0,z>0)
    pdfyz = zeros(1,N^2);       %2D pdf (wall where y>0,z>0)
    
    for n =1:nTup %look at each data point separately
        
        dat = Data(n,:);
        
        %determine i,j,k indices at which to compute pdf values
            minind = floor((dat - xo - h)./ delta) -1;
            maxind = ceil((dat -xo + h)./ delta) +1;
            minind(minind<1 | isnan(minind))=1;
            maxind(maxind>N | isnan(maxind))=N;
        
        ICoords = Coords(1,minind(1):maxind(1));
        if dat(1)~=zval(1)
            ICoords(ICoords == Coords(1,mergecoord(1)))=[];
        end
        
        if dim >1
            JCoords = Coords(2,minind(2):maxind(2));
        if dat(2)~=zval(2)
            JCoords(JCoords == Coords(2,mergecoord(2)))=[];
        end
            if dim == 3
                KCoords = Coords(3,minind(3):maxind(3));
        if dat(3)~=zval(3)
            KCoords(KCoords == Coords(3,mergecoord(3)))=[];
        end
            end
        end
        

        %use c mex functions to calculate pdf
        if dim==1
            pdfadd = mdKDE_1d(dat,minind,maxind,ICoords,N,h);        
            pdfadd(isnan(pdfadd))=0;
            pdfcenter = pdfcenter+pdfadd;
            
        elseif dim ==2
            
            if dat(1)~=zval(1) && dat(2)~=zval(2)
                pdfadd = mdKDE_2d(dat,minind,maxind,ICoords,JCoords,N,h);
            else
                pdfadd=0;
            end
            
            pdfadd(isnan(pdfadd))=0;
            pdfcenter = pdfcenter+pdfadd;
            %marginal pdfs: pxo(x) for x~=0,y=0 and pyo(y) for x=0,y~=0
            
           % if atom(1)==1 && dat(1)==0 && dat(2)~=0 %pyo: zero effect for x
           if dat(1)==zval(1) && dat(2)~=zval(2) %pyo: zero effect for x
                pdfyadd = mdKDE_1d(dat(2),minind(2),maxind(2),JCoords,N,h1D(2));
                pdfyadd(isnan(pdfyadd))=0;
                pdfy = pdfy + pdfyadd;
            end
            
           % if atom(2)==1 && dat(1)~=0 && dat(2)==0 %pxo: zero effect for y
           if dat(1)~=zval(1) && dat(2)==zval(2) %pxo: zero effect for y
                pdfxadd = mdKDE_1d(dat(1),minind(1),maxind(1),ICoords,N,h1D(1));
                pdfxadd(isnan(pdfxadd))=0;
                pdfx = pdfx + pdfxadd;
            end
            
        elseif dim ==3
            
            if dat(1)~=zval(1) && dat(2)~=zval(2) && dat(3)~=zval(3) %center points
                %make sure Icoords, JCoords, KCoords do no include
                %zero-coordinate in each dimension
                
                
                pdfadd= mdKDE_3d(dat,minind,maxind,ICoords,JCoords,KCoords,N,h);
                pdfadd(isnan(pdfadd))=0;
                pdfcenter = pdfcenter+pdfadd;
            elseif dat(1)==zval(1) && dat(2)~=zval(2) && dat(3)~=zval(3) %p(y,z)
                pdfyzadd = mdKDE_2d([dat(2) dat(3)],[minind(2) minind(3)],...
                    [maxind(2) maxind(3)],JCoords,KCoords,N,[h2D(2) h2D(3)]);
                pdfyzadd(isnan(pdfyzadd))=0;
                pdfyz = pdfyz+pdfyzadd;
            elseif dat(1)~=zval(1) && dat(2)==zval(2)  && dat(3)~=zval(3) %p(x,z)
                pdfxzadd = mdKDE_2d([dat(1) dat(3)],[minind(1) minind(3)],...
                    [maxind(1) maxind(3)],ICoords,KCoords,N,[h2D(1) h2D(3)]);
                pdfxzadd(isnan(pdfxzadd))=0;
                pdfxz = pdfxz+pdfxzadd;
            elseif dat(1)~=zval(1) && dat(2)~=zval(2)  && dat(3)==zval(3) %p(x,y)
                pdfxyadd = mdKDE_2d([dat(1) dat(2)],[minind(1) minind(2)],...
                    [maxind(1) maxind(2)],ICoords,JCoords,N,[h2D(1) h2D(2)]);
                pdfxyadd(isnan(pdfxyadd))=0;
                pdfxy = pdfxy+pdfxyadd;
            elseif dat(1)==zval(1) && dat(2)==zval(2) && dat(3)~=zval(3) %p(z)
                pdfzadd = mdKDE_1d(dat(3),minind(3),maxind(3),KCoords,N,h1D(3));
                pdfzadd(isnan(pdfzadd))=0;
                pdfz = pdfz + pdfzadd;
            elseif dat(1)==zval(1) && dat(2)~=zval(2) && dat(3)==zval(3) %p(y)
                pdfyadd = mdKDE_1d(dat(2),minind(2),maxind(2),JCoords,N,h1D(2));
                pdfyadd(isnan(pdfyadd))=0;
                pdfy = pdfy + pdfyadd;
            elseif dat(1)~=zval(1) && dat(2)==zval(2) && dat(3)==zval(3) %p(x)
                pdfxadd = mdKDE_1d(dat(1),minind(1),maxind(1),ICoords,N,h1D(1));
                pdfxadd(isnan(pdfxadd))=0;
                pdfx = pdfx + pdfxadd;
            end
            
        end
        
    end
      
    pdfcenter(isnan(pdfcenter))=0;

    
    if dim==1
        %pdfcenter = pdfcenter(2:end);
        pdfcenter(mergecoord(1))=0;
        pdfcenter = pdfcenter./sum(pdfcenter).*fx;
    
    elseif dim==2
        
        pdfcenter = reshape(pdfcenter,N,N); 
        %pdfcenter = pdfcenter(2:end,2:end);
        if sum(sum(pdfcenter))==0
            pdfcenter = ones(size(pdfcenter));
        end

        pdfcenter(mergecoord(1),:)=0;
        pdfcenter(:,mergecoord(2))=0;
        pdfcenter = (pdfcenter./(sum(sum(pdfcenter)))).*fxy;
        
        pdfy(isnan(pdfy))=0;
        %pdfy = pdfy(2:end);
        if sum(pdfy)==0
         pdfy = ones(size(pdfy));   
        end

        pdfy(mergecoord(2))=0;
        pdfy = (pdfy./sum(pdfy)).*fy;
        
        
        pdfx(isnan(pdfx))=0;
        %pdfx = pdfx(2:end);
        if sum(pdfx)==0
        pdfx=ones(size(pdfx));
        end
        pdfx(mergecoord(1))=0;
        pdfx = (pdfx./sum(pdfx)).*fx;
         
        pdfy(isnan(pdfy))=0;
        pdfx(isnan(pdfx))=0;
        pdfcenter(isnan(pdfcenter))=0;
               
    elseif dim==3
            
        pdfx(isnan(pdfx))=0;
        %pdfx = pdfx(2:end);
        if sum(pdfx)==0
        pdfx=ones(size(pdfx));
        end
        pdfx(mergecoord(1))=0;
        pdfx = (pdfx./sum(pdfx)).*fx;
        
        pdfy(isnan(pdfy))=0;
        %pdfy = pdfy(2:end);
        if sum(pdfy)==0
        pdfy=ones(size(pdfy));
        end
        pdfy(mergecoord(2))=0;
        pdfy = (pdfy./sum(pdfy)).*fy;
        
        pdfz(isnan(pdfz))=0;
        %pdfz = pdfz(2:end);
        if sum(pdfz)==0
        pdfz=ones(size(pdfz));
        end
        pdfz(mergecoord(3))=0;
        pdfz = (pdfz./sum(pdfz)).*fz;
        
        pdfxy(isnan(pdfxy))=0;
        pdfxy = reshape(pdfxy,N,N);
        %pdfxy = pdfxy(2:end,2:end);
        if sum(sum(pdfxy))==0
        pdfxy=ones(size(pdfxy));
        end
        pdfxy(mergecoord(1),:)=0;
        pdfxy(:,mergecoord(2))=0;
        pdfxy = (pdfxy./sum(sum((pdfxy)))).*fxy;
        
        pdfxz(isnan(pdfxz))=0;
        pdfxz = reshape(pdfxz,N,N);
        %pdfxz = pdfxz(2:end,2:end);
        if sum(sum(pdfxz))==0
        pdfxz=ones(size(pdfxz));
        end
        pdfxz(mergecoord(1),:)=0;
        pdfxz(:,mergecoord(3))=0;
        pdfxz = (pdfxz./sum(sum((pdfxz)))).*fxz;
        
        pdfyz(isnan(pdfyz))=0;
        pdfyz = reshape(pdfyz,N,N);
        %pdfyz = pdfyz(2:end,2:end);
        if sum(sum(pdfyz))==0
        pdfyz=ones(size(pdfyz));
        end
        pdfyz(mergecoord(2),:)=0;
        pdfyz(:,mergecoord(3))=0;
        pdfyz = (pdfyz./sum(sum((pdfyz)))).*fyz;
       
        pdfcenter = reshape(pdfcenter,N,N,N);
        %pdfcenter = pdfcenter(2:end,2:end,2:end);
        if sum(sum(sum(pdfcenter)))==0
        pdfcenter=ones(size(pdfcenter));
        end
        pdfcenter(mergecoord(1),:,:)=0;
        pdfcenter(:,mergecoord(2),:)=0;
        pdfcenter(:,:,mergecoord(3))=0;
        pdfcenter = pdfcenter./sum(sum(sum(pdfcenter))).*fxyz;
        
    end
    
    %merge the marginal and central pdfs together
    if dim==1 
        pdf = pdfcenter;
        pdf(mergecoord(1)) = f0;
    elseif dim==2 
        
        pdf = pdfcenter;
       % fprintf('step 1 center pdf = %6.4f\n',sum(sum(pdf)))

      %  fprintf('sum of y insert location = %6.4f, insert size = %6.4f\n',sum(pdf(mergecoord(1),:)),sum(pdfy))
        pdf(mergecoord(1),:)=pdfy';
      %  fprintf('step 2 center+y pdf = %6.4f\n',sum(sum(pdf)))
        
        
      %  fprintf('sum of x insert location = %6.4f, insert size = %6.4f\n',sum(pdf(:,mergecoord(2))),sum(pdfx))
        pdf(:,mergecoord(2))=pdfx;
      %  fprintf('step 3 center+x+y pdf = %6.4f\n',sum(sum(pdf)))
        
     %   fprintf('sum of 00 insert location = %6.4f, insert size = %6.4f\n',sum(pdf(mergecoord(1),mergecoord(2))),f0) 
        
        pdf(mergecoord(1),mergecoord(2))=f0;
     %   fprintf('step 4 center+y+x+middle pdf = %6.4f\n',sum(sum(pdf)))
        if sum(sum(pdf))<.999
        fprintf('sum of 2d pdf = %5.4f\n',sum(sum(pdf)))
        fprintf('y = %4.3f, x = %4.3f, zeros = %4.3f, center = %4.3f\n', sum(pdfy), sum(pdfx), f0,sum(sum(pdfcenter)))
        fprintf('fy = %5.3f, fx = %5.3f, f0=%5.3f, fxy=%5.3f\n',fy,fx,f0,fxy)
        
        fprintf('employing fixed binning instead\n')
        %use fixed bin pdf instead
        pdf = compute_pdf_fixedbins(Data,N,bin_scheme, Range);
        
        end
    elseif dim==3  
        pdf = pdfcenter;
        pdf(mergecoord(1),:,:)=pdfyz;
        pdf(:,mergecoord(2),:)=pdfxz;
        pdf(:,:,mergecoord(3))=pdfxy;
        pdf(mergecoord(1),mergecoord(2),:)=pdfz;
        pdf(mergecoord(1),:,mergecoord(3))=pdfy;
        pdf(:,mergecoord(2),mergecoord(3))=pdfx;
        pdf(mergecoord(1),mergecoord(2),mergecoord(3))=f0;

        if sum(sum(sum(pdf)))<.99
        sprintf('sum of 3d pdf = %5.4f\n',sum(sum(sum(pdf))))
        sprintf('h1 = %4.3f, h2 = %4.3f, h3=%4.3f\n',h(1),h(2),h(3))
         fprintf('employing fixed binning instead\n')
        end
        
        %use fixed bin pdf instead
        pdf = compute_pdf_fixedbins(Data,N,bin_scheme, Range);
        
    end
    
elseif strcmp(method,'fixed')==1 %%%%%%%%%%% Fixed Bin Method  %%%%%%%%%%%%%%%%%
    
 
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

