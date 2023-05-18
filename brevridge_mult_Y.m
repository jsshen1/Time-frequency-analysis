function [Cs] = brevridge_mult_Y(Tx, fre, delta, clwin)
% Function brevridge_mult : extracts the ridges of a multicomponent signal
% Inputs:
%   Tx TF transform of s (SST, STFT, RM...)
%   lfs : log(fs) frequencies
%   nr : number of ridges
%   lambda : lambda parameter
%   clwin : frequency clearing window
% Outputs:
%   Cs Cs(:,j) is the j-th ridge location

if nargin<2
    fre=1:2;
    delta=0.8;
    clwin=5;
elseif nargin<3
    delta=0.8;
    %lambda=0.001;
    clwin=5;
elseif nargin<4
    clwin=5;
end

Txs = Tx;

[na,N] = size(Txs);

%Cs = zeros(N, nr);
%Es = zeros(nr, 1);
flag=1;
k=1;
%Cs=[];
max_e=N;
%max_e=0;
while flag==1
   [Cs_temp] = brevridge_Y(Txs, fre, 0.001);
   
   e=0;
for i=1:N
e=e+Txs(Cs_temp(i),i);
end

Cs(:,k)=Cs_temp;

if e>max_e
max_e=e;
end

%if e<N*delta;
if e<max_e*delta;
    Cs(:,k)=[];
flag=0;
break;
end

 

  for b=1:N
        Txs(max(1,Cs_temp(b)-clwin):min(na,Cs_temp(b)+clwin),b)=0;
  end
 k=k+1;   
  
end

% 
% for j=1:nr
%     [Cs(:,j), Es(j)] = brevridge(Txs, fs, lambda);
% 
%     % Remove this curve from the representation
%     % Max/min frequencies for each time step
%     for b=1:N
%         Txs(max(1,Cs(b,j)-clwin):min(na,Cs(b,j)+clwin),b)=0;
%     end
% end

Cs = Cs';

end
