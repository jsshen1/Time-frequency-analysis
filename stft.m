function tfr = stft(x,hlength,sigma)
% Computes the SST (Ts)  of the signal x.
% INPUT
%    x      :  Signal needed to be column vector.
%    hlength:  The hlength of window function.
% OUTPUT
%    Ts     :  The SST

[xrow,xcol] = size(x);

if (xcol~=1),
 error('X must be column vector');
end; 

if (nargin < 1),
error('At least 1 parameter is required');
end;

if (nargin < 2),
hlength=round(xrow/5);
le=round(hlength/2);
end;

if (nargin < 3),
le=round(hlength/2);
end;

%Siglength=xrow;
hlength=hlength+1-rem(hlength,2);
ht = linspace(-0.5,0.5,hlength);
ht=ht';


% Gaussian window
h = exp(-pi/sigma^2*ht.^2);
[hrow,~]=size(h); 
Lh=(hrow-1)/2; 

N=xrow;
t=1:xrow;
[~,tcol] = size(t);

tfr1= zeros (N,tcol) ; 

for icol=1:tcol
ti= t(icol); 
tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);%´°µÄ·¶Î§
indices= rem(N+tau,N)+1;
rSig = x(ti+tau,1);
%rSig = hilbert(real(rSig));
tfr1(indices,icol)=rSig.*conj(h(Lh+1+tau));
%tfr2(indices,icol)=rSig.*conj(dh(Lh+1+tau));
end

tfr1=fft(tfr1);
%tfr2=fft(tfr2);

tfr=tfr1(1:round(N/2),:);

end