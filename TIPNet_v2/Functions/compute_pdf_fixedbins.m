function [pdf, Coords]= compute_pdf_fixedbins(Data,N,bin_scheme, Range)
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

%measures to determine at which grid points to consider contributions from
%each data point: delta, h, xo

Coords = zeros(dim,N);
Edges = zeros(dim,N+1);

%compute non-zero fraction for each dimension to determine whether to
%correct for zero-effect (if more than 10% zero values)

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




