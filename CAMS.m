function [xn,peaks,CoCe,shiftvalue] = CAMS(reftic,refmz,xtic,xmz,phi,max_shift,valve)
%  CAMS	
%  Input 
%         reftic: reference chromatogram (row vector)
%         refmz: matrix with mass spectra of reference chromatogram data
%         xtic: chromatogram to be aligned (row vector)
%         xmz:  matrix with mass spectra of chromatogram to be aligned
%         phi:  level for DWT
%         max_shift: allowed maximum shift parameter for each segment
%         valve:the valve of cross correlation of mass spectrums
%  Output
%         xn: aligned chromatogram
%         peaks: detected peaks
%         CoCe: corrcoef correlation of two chromatograms
%         shiftvalue: shift value of each peak
%  Examples:
%         please see testCAMS.m
%
%  Main reference:
%         (1) Wong, Jason W. H., Durante, Caterina, and Cartwright, Hugh M., Application of Fast Fourier Transform Cross-Correlation for the Alignment of Large Chromatographic and Spectral Datasets. Analytical Chemistry 77 (17), 5655 (2005).
%         (2) Veselkov, K. A. et al., Recursive Segment-Wise Peak Alignment of Biological H-1 NMR Spectra for Improved Metabolic Biomarker Recovery. Analytical Chemistry 81 (1), 56 (2009).
%         (3) Zhang, Z. M., Chen, S., and Liang, Y. Z., Peak alignment using wavelet pattern matching and differential evolution. Talanta 83 (4), 1108 (2011).
%         (4) Nielsen, Niels-Peter Vest, Carstensen, Jens Michael, and Smedsgaard, J?rn, Aligning of single and multiple wavelength chromatographic profiles for chemometric data analysis using correlation optimised warping. Journal of Chromatography A 805 (1-2), 17 (1998).
%
%  Yi-Bao Zheng @ central south university on April 25,2012
%  E-mail: matchcoder@gmail.com
%  Web: www.chemosolv.com/blog/matchcoder
%
if nargin<7; valve=0.8; end;
if nargin<6; max_shift=90; end;
if nargin<5; phi=6; end;
if (nargin<4)
   help align;
   return;
end
l=6;
[xs,xln]=wavelet_denoising(xtic,l,6,'soft');
[ss,sln]=wavelet_denoising(reftic,l,6,'soft');
peaks=peak_detection(xtic,3,xln,phi);
peaks=peak_clustering(xtic,peaks,1);
xs=xtic;
xn=xs;
shiftvalue=zeros(1,size(peaks,1));
for i=size(peaks,1):-1:1
    si=max(1,peaks(i,1)-max_shift);
    ei=min(length(xs),peaks(i,3)+max_shift);
    [lag,fccc]=fftcc(xs(si:ei),ss(si:ei),max_shift); 
    
    for j=length(lag):-1:1
        if i>2 && i<size(peaks,1)
            if (lag(j)>0) && (lag(j)>(peaks(i+1,2)-peaks(i,2)))
                lag(j) =[];
            end
            if (lag(j)<0) && (abs(lag(j))> (peaks(i,2)-peaks(i-1,2)))
                lag(j) =[];
            end     
        end
    end
    for j=1:length(lag)
        if (peaks(i,2)+lag(j))>0 && (peaks(i,2)+lag(j)<length(reftic))
            R(j) = corr2(xmz(:,(peaks(i,2)-1):(peaks(i,2)+1)),refmz(:,(peaks(i,2)+lag(j)-1):(peaks(i,2)+lag(j)+1)));
        end
    end
    if max(R)>=valve
        [value,shift] = max(R);
        shiftvalue(i)=lag(shift);
    else
       peaks(i,:) = [];
       shiftvalue(i) = [];
    end
    CoCe(i)=max(R);
    R=[];
end

for i=1:size(peaks,1)  
    if i==1
        if peaks(i,1)+shiftvalue(i) <= 1
            xn(1:(peaks(i,3)+shiftvalue(i))) = xs(abs(shiftvalue(i)-1):peaks(i,3));
        elseif peaks(i,1)==1
            xn(1:peaks(i,3)) = xs(1:peaks(i,3));
        else
            xt=interp1(1:peaks(i,1), xs(1:peaks(i,1)),linspace(1,peaks(i,1),peaks(i,1)+shiftvalue(i)));
            xn(1:(peaks(i,1)+shiftvalue(i))) = xt;
            xn((peaks(i,1):peaks(i,3))+shiftvalue(i)) = xs(peaks(i,1):peaks(i,3)); 
        end
        if size(peaks,1)==1
            xt=interp1(peaks(i,3):length(xs), xs(peaks(i,3):length(xs)),linspace(peaks(i,3),length(xs),length(xs)-peaks(i,3)-shiftvalue(i)+1));
            xn((peaks(i,3)+shiftvalue(i)):length(xs)) = xt;
        end
    elseif i==size(peaks,1)
        if peaks(i,3)+shiftvalue(i) < length(xn)
            xt=interp1(peaks(i,3):length(xs), xs(peaks(i,3):length(xs)),linspace(peaks(i,3),length(xs),length(xs)-peaks(i,3)-shiftvalue(i)+1));
            xn((peaks(i,3)+shiftvalue(i)):length(xs)) = xt;
            xt=interp1(peaks(i-1,3):peaks(i,1), xs(peaks(i-1,3):peaks(i,1)),linspace(peaks(i-1,3),peaks(i,1),peaks(i,1)+shiftvalue(i)+1-(peaks(i-1,3)+shiftvalue(i-1))));
            xn((peaks(i-1,3)+shiftvalue(i-1)):(peaks(i,1)+shiftvalue(i))) = xt;
            xn((peaks(i,1):peaks(i,3))+shiftvalue(i)) = xs(peaks(i,1):peaks(i,3));
        end
    else
        if peaks(i-1,3)+shiftvalue(i-1)>peaks(i,1)+shiftvalue(i)
             xn((peaks(i,1)+shiftvalue(i)):(peaks(i-1,3)+shiftvalue(i-1))) =  xn((peaks(i,1)+shiftvalue(i)):(peaks(i-1,3)+shiftvalue(i-1))) + xs(peaks(i,1):(peaks(i-1,3)+shiftvalue(i-1)-shiftvalue(i)));
             xn((peaks(i-1,3)+shiftvalue(i-1)):(peaks(i,3)+shiftvalue(i))) =  xs((peaks(i-1,3)+shiftvalue(i-1)-shiftvalue(i)):peaks(i,3));
        else
            xt=interp1(peaks(i-1,3):peaks(i,1), xs(peaks(i-1,3):peaks(i,1)),linspace(peaks(i-1,3),peaks(i,1),peaks(i,1)+shiftvalue(i)+1-(peaks(i-1,3)+shiftvalue(i-1))));
            xn((peaks(i-1,3)+shiftvalue(i-1)):(peaks(i,1)+shiftvalue(i))) = xt;
            xn((peaks(i,1):peaks(i,3))+shiftvalue(i)) = xs(peaks(i,1):peaks(i,3));
       end
    end
end