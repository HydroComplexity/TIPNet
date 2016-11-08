function [ry, py, plotrange, ffilter]=...
    Filter_function(xvector,yvector,centerfrequency,filterwidth,filtershape,FILTERMODE)
% Separate graph windows for the original and filtered signals.
% Computes and plots fourier filter for signal yvector.  
% Centerfrequency and filterwidth are the center frequency and
% width of the pass band, in harmonics, 'filtershape' determines
% the sharpness of the cut-off. Plot modes: mode=1 linear x and y; 
% mode=2 log x linear y; mode=3 linear x; log y; mode=4 log y log x

fy=fft(yvector);
lft1=[1:(length(fy)/2)];
lft2=[(length(fy)/2+1):length(fy)];
% Compute filter shape.
if strcmp(FILTERMODE,'Band-pass'),
       ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
       ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
       ffilter=[ffilter1,ffilter2];
end
if strcmp(FILTERMODE,'High-pass')
       centerfrequency=length(xvector)/2;
       ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
       ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
       ffilter=[ffilter1,ffilter2];
end
if strcmp(FILTERMODE,'Low-pass')
       centerfrequency=0;
       ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
       ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
       ffilter=[ffilter1,ffilter2];
end
if strcmp(FILTERMODE,'Band-reject (notch)')
       ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
       ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
       ffilter=1-[ffilter1,ffilter2]; 
end
if strcmp(FILTERMODE,'Comb pass')
    n=2;
    ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
    ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
    while n<30,
       ffilter1=ffilter1+shape(lft1,n*(centerfrequency+1),filterwidth,filtershape);
       ffilter2=ffilter2+shape(lft2,length(fy)-n*(centerfrequency+1),filterwidth,filtershape);
       n=n+1;
    end
       ffilter=[ffilter1,ffilter2];
end
if strcmp(FILTERMODE,'Comb notch')
    n=2;
    ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
    ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
    while n<30,
       ffilter1=ffilter1+shape(lft1,n*(centerfrequency+1),filterwidth,filtershape);
       ffilter2=ffilter2+shape(lft2,length(fy)-n*(centerfrequency+1),filterwidth,filtershape);
       n=n+1;
    end
       ffilter=1-[ffilter1,ffilter2];
end

if length(fy)>length(ffilter), ffilter=[ffilter ffilter(1)];end
ffy=fy.*ffilter;  % Multiply filter by Fourier Transform of signal
ry=real(ifft(ffy));

py=fy .* conj(fy); % Compute power spectrum
plotrange=2:length(fy)/2;


function g = shape(x,pos,wid,n)
%  shape(x,pos,wid,n) = peak centered on x=pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  Shape is Lorentzian (1/x^2) when n=0, Gaussian (exp(-x^2))
%  when n=1, and becomes more rectangular as n increases.
%  Example: shape([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
if n==0
    g=ones(size(x))./(1+((x-pos)./(0.5.*wid)).^2);
else
    g = exp(-((x-pos)./(0.6.*wid)) .^(2*round(n)));
end

function r = range(arr)
r = max(arr) - min(arr);