function [pos,dx]=cwt_local_maximum(x,hpw)

dx = -cwt(x,2*hpw,'haar');
pos=find(sign(-sign(diff(sign(dx))+0.5)+1))';

for i=1:size(pos,1)
    if x(pos(i)+1)>x(pos(i))
       pos(i)=pos(i)+1;
    end
end

% for i =1:hpw
%     if x(i)>x(pos(1))
%         pos(i) = i;
%     end  
% end
% 
% for i = size(x,2):-1:(size(x,2)-hpw+1)
%     if x(i)>x(pos(1))
%         pos(i) = i;
%     end  
% end

