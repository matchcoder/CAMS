%load the data
load Data247
load Data260

% baseline correction using airPLS
lambda = 1e4;
[xbc,xb]=airPLS(Data247.TIC, lambda,2,0.05);
[sbc,sb]=airPLS(Data260.TIC, lambda,2,0.05);

%alignment
tic
[xn,peaks,CoCe,shiftvalue] = CAMS(sbc,Data260.MZ,xbc,Data247.MZ,6,20,0.8);
toc

%corrcoef correlation of two chromatograms before alignment
corrcoef(xbc,sbc)
%corrcoef correlation of two chromatograms after alignment
corrcoef(xn,sbc)

%figure out the result
figure
plot(sbc+500000,'b')
hold on
plot(xn+250000,'r')
hold on
plot(xbc,'g')
legend('Reference','After aligned','To be aligned')
