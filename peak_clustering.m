function peaks_c=peak_clustering(x,peaks,n)

nPeakNum=size(peaks,1);
vGaps=zeros(nPeakNum-1,1);
peaks_c=peaks;

for i=1:(nPeakNum-1)
    vGaps(i)=peaks_c(i+1,1)-peaks_c(i,3);
    if vGaps(i)<=n
        peaks_c(i+1,1)=peaks_c(i,1);
        if x(peaks_c(i+1,2)) < x(peaks_c(i,2))
            peaks_c(i+1,2) = peaks_c(i,2);
        end
    end
end

vSmallGapsInd=find(vGaps<=n);

for i=1:size(vSmallGapsInd,1)
    peaks_c(vSmallGapsInd(i)-i+1,:)=[];
end

