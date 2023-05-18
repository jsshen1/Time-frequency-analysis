function [extr_Sig,ampmatrix,iniIFfit] = ICCD(Sig,SampFreq,iniIFset,orderIF,orderamp,alpha)
% Intrinsic Chirp Component Decomposition(ICCD)
% the code is only for complex-valued data analysis
%%%%%%%%%%%%%%%%%%%%%%%  input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sig£ºmeasured signal,a row vector
% SampFreq: sampling frequency
% iniIFset: the instantaneous frequency series
% orderIF: the order of the Fourier model used for fitting the instantaneous frequency of the signal
% orderamp£ºFourier order for characterizing signal amplitudes
% alpha£ºTikhonov regularization parameter for ICCD.
%%%%%%%%%%%%%%%%%%%%%%%  output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extr_Sig: the reconstructed components
% ampmatrix: the estimated amplitudes
% iniIFfit£ºthe fitted instantaneous frequencies (IFs)

if (isreal(Sig))
Sig = hilbert(Sig);
end

[multin,N] = size(iniIFset);%multin denotes the number of the components£»N is the number of the samples
% when multin is larger than one, the algorithm performs the joint-optimization scheme.


dt = [0:N-1]/SampFreq; %time

phase = zeros(multin,N);
iniIFfit = zeros(multin,N);
for i = 1:multin
 [temp1,temp2] = coef_ovefour(iniIFset(i,:),SampFreq,orderIF);  % fitting the input IFs with a Fourier model
 iniIFfit(i,:) = temp1;%fitted IF
 phase(i,:) = temp2; % fitted phase 
end
%%%%%%%%%%%%%%%%%%%%Construction of the Fourier matrix%%%%%%%%%%%%%%%%%%%%%%%
f0 = SampFreq/2/N;%%Base Frequency
l_t = 2*orderamp + 1;%%the number of the Fourier coefficients
tmatrix = zeros(N,l_t);
tmatrix(:,1) = ones(N,1);
for j = 2:l_t
        tmatrix(:,j) = cos(2*pi*f0*(j-1)*dt); 
        if j >(l_t+1)/2
            tmatrix(:,j) = sin(2*pi*f0*(j-((l_t+1)/2))*dt);
        end
end

l_kernel = multin*l_t;%
kmatrix = zeros(N,l_kernel);%
for i = 1:multin
   C = exp(sqrt(-1)*2*pi*phase(i,:));%
   kmatrix(:,((i-1)*l_t+1):i*l_t) =  spdiags(C(:), 0, N, N)*tmatrix;
end

Sig = Sig(:);
%% regularized least-squares solution
Imatrix = speye(size(kmatrix,2));%Identity matrix
theta =(alpha*Imatrix + kmatrix'*kmatrix)\(kmatrix'*Sig);%coefficient vector
%% extract each component
extr_Sig = zeros(multin,N);
for i = 1 : multin 
    extr_Sig(i,:) = kmatrix(:,(i-1)*l_t+1:i*l_t)*theta((i-1)*l_t+1:i*l_t);
end
%% extract amplitude %%%%
complex_amp = zeros(multin,N);
for i =  1 :multin
    complex_amp(i,:) = tmatrix*theta(((i-1)*l_t+1):(i*l_t));
end
ampmatrix = abs(complex_amp);
end
    


