function [xs,xln]=wavelet_denoising(x,l,t,type)

nx=size(x,2);
[xW,xL] = wavedec(x,l,'db8');
T = t*median(abs(xW(floor(nx/2)+1:nx)))/.6745;
xWT = thresholding(xW,T,type);
xs = waverec(xWT,xL,'db8');
xn=x-xs;
xln=local_noise_estimate(xn,41);
%xs=abs(xs);

function y = thresholding(x, t, type)

switch lower(type)
    case 'hard'
        y = hard_thresholding(x,t);
    case 'soft'
        y = soft_thresholding(x,t);
    case 'sure'
        y = sure_thresholding(x,t);
    otherwise
        error('Unkwnown.');
end

function y = hard_thresholding(x,t)

y = (abs(x) > t).*x; 

function y = sure_thresholding(x,t)

y = x .* max( 1-t^2 ./ abs(x).^2, 0 );

function y = soft_thresholding(x,t)

y = sign(x).*(abs(x) >= t).*(abs(x) - t); 

function xln=local_noise_estimate(xn,wsize)

halfwsize = (wsize - 1) / 2;
xnp = [xn(halfwsize:-1:1) xn xn(end:-1:end-halfwsize)];
for i=1:length(xn)
   xln(i) = mad(xnp(i:i+wsize)) / .67;
end
