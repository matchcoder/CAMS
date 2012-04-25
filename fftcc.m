function [lag,fccc]= fftcc(x,s,shift,hpw)
if nargin < 4
    hpw=3;
end
n=size(s,2);
shift=max(1,min(n,shift));
S=fft(s);
X=fft(x);
fccc=real(ifft(S.*conj(X)./n));
lag=cwt_local_maximum(fccc,hpw);
%lag=lag(((lag>hpw)&(lag<=shift))|((lag>=n-shift)&(lag<n-hpw)));

[max_l,ind_l]=max(fccc(1:hpw+1));
[max_r,ind_r]=max(fliplr(fccc(n-hpw-1:n)));
%lag=[lag;(1:hpw)';ind_l;(n-ind_r)];
lag=[lag;ind_l;(n-ind_r)];
lag=unique(lag);
%sfccc=fccc(lag);
lag=lag-1;
if any(lag>n/2)
   lag(lag>n/2)=lag(lag>n/2)-n; 
end