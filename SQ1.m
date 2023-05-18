function [tfr1,SST,Dtfr,omega,T] = SQ1(t,sig,fre,wlen,sigma)
%
tLen = length(t);
N = length(fre);%频率总数量,freqlow=0，alpha=1，freqhigh=255
sig=sig(:);
alpha=fre(2)-fre(1);
tfr = stft(sig,wlen,sigma);
tfr1=tfr;
tfr=abs(tfr);
Dtfr=[tfr(2:end,:)-tfr(1:end-1,:);tfr(end,:)-tfr(end-1,:)]/(alpha);
Dtfr1=[Dtfr(2:end,:); Dtfr(end,:)];
T=Dtfr.*Dtfr1;

omega=zeros(N,tLen);
for i = 1:tLen
    loc=find(T(:,i)<0);
    omega(loc,i)=fre(loc);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    	%% Synchro-squeezing transform
% tfr=abs(tfr);
SST = zeros(N,tLen); 
Ex = mean(tfr.^2);
Threshold = 0.00001*Ex;	% originally it was 1e-6*Ex
% Threshold = 0;
for i = 1:tLen
    for j = 1: N
        if abs(tfr(j,i)) > Threshold
            k =1+round((omega(j,i)-fre(1))/alpha);
            if k>=1 && k<=N% Vertical reassignment: SST 
                SST(k,i) = SST(k,i) +tfr1(j,i);
            end
        end
    end
end

end
