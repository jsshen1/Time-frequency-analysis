function [fitf,finte] = coef_ovefour(f,SampFreq,orderamp)
%smooth the discrete series f by least-squares fitting with a Fourier model£ºsee paper£º
% Chen S, Peng Z, Yang Y, et al. Intrinsic chirp component decomposition by using Fourier Series representation[J]. Signal Processing, 2017, 137: 319-327.
% f£ºmeasured discrete series
% SampFreq£ºsampling frequency
% orderamp£ºFourier order of the Fourier model
% fitf£ºfitting result
% finte£ºthe integral of fitf
f = f(:);
N = length(f);
dt = [0:N-1]/SampFreq; %time
%%%%%%%%%%%%%%%%%%%%Construction of the Fourier matrix%%%%%%%%%%%%%%%%%%%%%%%
f0 = SampFreq/2/N;%Base Frequency
orderamp = 2*orderamp + 1;%the number of the Fourier coefficients
tmatrix = zeros(N,orderamp);
tmatrix(:,1) = ones(N,1);%
for j = 2:orderamp
        tmatrix(:,j) = cos(2*pi*f0*(j-1)*dt); 
        if j >(orderamp+1)/2
            
            tmatrix(:,j) = sin(2*pi*f0*(j-((orderamp+1)/2))*dt);
        end
end
%%%%%%%%%%%%%%%%%%%%%%%%The integral of the Fourier model above%%%%%%%%%%%%%%%%%%%%%%
tmatrix_inte = zeros(N,orderamp);
tmatrix_inte(:,1) = dt;
for j = 2:orderamp
        tmatrix_inte(:,j) = 1/(2*pi*f0*(j-1))*sin(2*pi*f0*(j-1)*dt); 
        if j >(orderamp+1)/2
            
            tmatrix_inte(:,j) = -1/(2*pi*f0*(j-((orderamp+1)/2)))*cos(2*pi*f0*(j-((orderamp+1)/2))*dt);
        end
end
alpha = 0.00005;%Tikhonov regularization parameter
Imatrix = eye(size(tmatrix,2));
coeff =(alpha*Imatrix + tmatrix'*tmatrix)\(tmatrix'*f);

fitf = tmatrix*coeff;
finte = tmatrix_inte*coeff;







