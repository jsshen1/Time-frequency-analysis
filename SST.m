function [Ts,omega,IF] = SST(x,hlength);
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
end;

Siglength=xrow;
hlength=hlength+1-rem(hlength,2);
ht = linspace(-0.5,0.5,hlength);ht=ht';

% Gaussian window
h = exp(-pi/0.32^2*ht.^2);
% derivative of window
dh = -2*pi/0.32^2*ht .* h; % g'

[hrow,hcol]=size(h); Lh=(hrow-1)/2; 

N=xrow;
t=1:xrow;

[trow,tcol] = size(t);


tfr1= zeros (N,tcol) ; 
tfr2= zeros (N,tcol) ; 

tfr= zeros (round(N/2),tcol) ; 
Ts= zeros (round(N/2),tcol) ; 

for icol=1:tcol,
ti= t(icol); tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
indices= rem(N+tau,N)+1; 
rSig = x(ti+tau,1);
%rSig = hilbert(real(rSig));
tfr1(indices,icol)=rSig.*conj(h(Lh+1+tau));
tfr2(indices,icol)=rSig.*conj(dh(Lh+1+tau));
end;

tfr1=fft(tfr1);
tfr2=fft(tfr2);

tfr1=tfr1(1:round(N/2),:);
tfr2=tfr2(1:round(N/2),:);

ft = 1:round(N/2);
bt = 1:N;

%%operator omega
nb = length(bt);
neta = length(ft);

va=N/hlength;
omega = zeros (round(N/2),tcol);
IF = zeros (round(N/2),tcol);
for b=1:nb
omega(:,b) = (ft-1)'+real(va*1i*tfr2(ft,b)/2/pi./tfr1(ft,b));
end 


omega=round(omega);

for b=1:nb%time
    % Reassignment step
    for eta=1:neta%frequency
        if abs(tfr1(eta,b))>0.0001%you can set much lower value than this.
            k = omega(eta,b);
            if k>=1 && k<=neta
                Ts(k,b) = Ts(k,b) + tfr1(eta,b);
            end
            if abs(omega(eta,b)-eta)<0.5%default frequency resolution is 1Hz.
                IF(eta,b)=1;
            end
        end
    end
end
Ts=Ts/(xrow/2);
end