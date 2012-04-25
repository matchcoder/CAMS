function peaks=peak_detection(x,hpw,xln,phi)
% x is the vector of chromatogram intensities.
% hpw is the minimum half peak width
% peak_points is the vector of the peak points detected in "x"

%x=x+(phi/3)*sign(rand(size(xln,1),size(xln,2))-0.5).*xln;

% peak is a local maximum of N neighboring points.
[pos,dx]=peak_position(x,hpw);

% elimation of the peaks with intensity lower than the noise level
SNR=x./xln;
[SNR,xb]=airPLS(SNR,5e2,2,0.05);
for i=length(pos):-1:1
    if abs(SNR(pos(i))) < phi
        pos(i)=[];
    end
end

% estimation peaks width
pos_all=find(sign(diff(sign(dx))));

ind=find(pos_all==pos(1));
if ind==1
    widths(i,:)=[1 pos_all(ind+1)];    
else
    widths(i,:)=[pos_all(ind-1) pos_all(ind+1)];
end
for i=2:1:length(pos)-1
    ind=find(pos_all==pos(i));
    widths(i,:)=[pos_all(ind-1) pos_all(ind+1)];
end
ind=find(pos_all==pos(end));
if ind==length(pos_all)
    widths(length(pos),:)=[pos_all(ind-1) length(x)];    
else
    widths(length(pos),:)=[pos_all(ind-1) pos_all(ind+1)];
end

for i=1:size(pos,1)
    if x(pos(i)+1)>x(pos(i))
       pos(i)=pos(i)+1;
    end
end

peaks=[widths(:,1) pos widths(:,2)];

for i=size(peaks,1):-1:2
    vGaps=peaks(i,1)-peaks(i-1,3);
    if vGaps == 0 
        if (x(peaks(i,2))-x(peaks(i-1,2)))< max(x(peaks(i,2)),x(peaks(i-1,2)))/10 
            if x(peaks(i,2))-x(peaks(i-1,2)) < 0
                 peaks(i,:) = [];
            elseif x(peaks(i,2))-x(peaks(i-1,2)) > 0
                 peaks(i-1,:) = [];
            end
        end
    end    
end


    



