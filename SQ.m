function SST = SQ(fre,tfr,IF)
%
[N,tLen]=size(tfr);
omega=zeros(N,tLen);
df=fre(2)-fre(1);
[n,~]=size(IF);
for i = 1:tLen
    for j=1:n
        loc=round(IF(j,i));
        omega(loc,i)=fre(loc);
    end
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    	%% Synchro-squeezing transform
% tfr=abs(tfr);
SST = zeros(N,tLen); 
Threshold=1e-6;
for i = 1:tLen
    for j = 1: N
        if abs(tfr(j,i)) > Threshold
            k =1+round((omega(j,i)-fre(1))/df);
            if k>=1 && k<=N% Vertical reassignment: SST 
                SST(k,i) = SST(k,i) +tfr(j,i);
            end
        end
    end
end

end
