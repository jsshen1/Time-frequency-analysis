function [SST,SST1,Dtfr,omega,IF,T] = SQ2(t,tfr1,fre,gamma)
%
tLen = length(t);
N = length(fre);%Ƶ��������,freqlow=0��alpha=1��freqhigh=255
alpha=fre(2)-fre(1);
tfr=abs(tfr1);%�˴�ΪSTFTģ
gamma1=gamma*max(max(tfr));
tfr((abs(tfr) < gamma1)) = 0;
Dtfr=[tfr(2:end,:)-tfr(1:end-1,:);tfr(end,:)-tfr(end-1,:)]/(alpha);
Dtfr1=[Dtfr(2:end,:); Dtfr(end,:)];

T=Dtfr.*Dtfr1;
omega=zeros(N,tLen);
IF=zeros(N,tLen);
for i = 1:tLen
    loc=find(T(:,i)<0);
    omega(loc,i)=fre(loc);
    IF(loc,i)=ones(length(loc),1);
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    	%%squeezing

SST = zeros(N,tLen); 
SST1 = zeros(N,tLen);
Ex = mean(tfr.^2);
Threshold = 0.0001*Ex;	%̫С�Ļ��������ܵ�����
for i = 1:tLen
    for j = 1:N
        if abs(tfr1(j,i)) > Threshold
            k =1+round((omega(j,i)-fre(1))/alpha);
            if k>=1 && k<=N% Vertical reassignment: SST 
                %����fre(k)����Ϊ����ÿ�������������
                SST(k,i) = SST(k,i) +tfr(j,i)*(fre(end)-fre(k)+1)*fre(k);
                SST1(k,i) = SST1(k,i)+tfr1(j,i);
            end
        end
    end
end
SST=SST/(tLen/2);
SST1=SST1/(tLen/2);
end
