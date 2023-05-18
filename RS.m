function [tfr,IF] = RS(x,hlength);
%   Reassignment transform
%	x       : Signal.
%	hlength : Window length.

%	tfr   : Time-Frequency Representation.

%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

[xrow,xcol] = size(x);

N=xrow;

if (xcol~=1),
 error('X must be column vector');
end;

if (nargin < 2),
 hlength=round(xrow/8);
end;

t=1:N;
ft = 1:round(N/2);

[trow,tcol] = size(t);

hlength=hlength+1-rem(hlength,2);
ht = linspace(-0.5,0.5,hlength);ht=ht';

% Gaussian window
h = exp(-pi/0.32^2*ht.^2);
% derivative of window
dh = -2*pi/0.32^2*ht .* h; % g'
%
th=h.*ht;

[hrow,hcol]=size(h); Lh=(hrow-1)/2;

tfr1= zeros (N,tcol);
tfr2= zeros (N,tcol);
tfr3= zeros (N,tcol);

va=N/hlength;
    
for icol=1:tcol,
ti= t(icol); tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
indices= rem(N+tau,N)+1;
rSig = x(ti+tau,1);
tfr1(indices,icol)=rSig.*conj(h(Lh+1+tau));
tfr2(indices,icol)=rSig.*conj(dh(Lh+1+tau));
tfr3(indices,icol)=rSig.*conj(th(Lh+1+tau));
end;

tfr1=fft(tfr1);
tfr2=fft(tfr2);
tfr3=fft(tfr3);

tfr1=tfr1(1:round(N/2),:);
tfr2=tfr2(1:round(N/2),:);
tfr3=tfr3(1:round(N/2),:);

omega = zeros(round(N/2),tcol);
omega2= zeros(round(N/2),tcol);

for b=1:N
omega1(:,b) = (ft-1)'+real(va*1i*tfr2(ft,b)/2/pi./tfr1(ft,b));
end

for a=1:round(N/2)
omega2(a,:) = t+(hlength-1)*real(tfr3(a,t)./tfr1(a,t));
end

omega1=round(omega1);
omega2=round(omega2);

Ts = zeros(round(N/2),tcol);
IF = zeros(round(N/2),tcol);
% Reassignment step
for b=1:N%time
    for eta=1:round(N/2)%frequency
        if abs(tfr1(eta,b))>0.0001
            k1 = omega1(eta,b);
            k2 = omega2(eta,b);
            if k1>=1 && k1<=round(N/2) && k2>=1 && k2<=N
                Ts(k1,k2) = Ts(k1,k2) + abs(tfr1(eta,b));
            end
            if abs(omega1(eta,b)-eta)<0.5%default frequency resolution is 1Hz.
                IF(eta,b)=1;
            end
        end
    end
end

tfr=Ts/(xrow/2);
