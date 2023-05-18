function [Te,IF, tfr] = SET(x,hlength)
%   Synchroextracting Transform
%	x       : Signal.
%	hlength : Window length.

%   IF   : Synchroextracting operator representation.
%   Te   : SET result.
%   tfr  : STFT result
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
%
%   Written by YuGang in Shandong University at 2016.5.13.

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

[hrow,hcol]=size(h); Lh=(hrow-1)/2;

tfr1= zeros (N,tcol);
tfr2= zeros (N,tcol);

for icol=1:tcol,
ti= t(icol); tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
indices= rem(N+tau,N)+1;
rSig = x(ti+tau,1);
tfr1(indices,icol)=rSig.*conj(h(Lh+1+tau));
tfr2(indices,icol)=rSig.*conj(dh(Lh+1+tau));
end;

tfr1=fft(tfr1);
tfr2=fft(tfr2);

tfr1=tfr1(1:round(N/2),:);
tfr2=tfr2(1:round(N/2),:);

va=N/hlength;
IF=zeros(round(N/2),tcol);
%tfr=zeros(round(N/2),tcol);
E=mean(abs(x));

for i=1:round(N/2)%frequency
for j=1:N%time
     if abs(tfr1(i,j))>0.8*E%if you are interested in weak signals, you can delete this line.
         %if abs(1-real(va*1i*tfr2(i,j)/2/pi./tfr1(i,j)))<0.5
         if abs(-real(va*1i*tfr2(i,j)/2/pi./tfr1(i,j)))<0.5
         IF(i,j)=1;
         end
     end
end
end
tfr=tfr1/(sum(h)/2);%the amplitude of tfr result has been pre-rectified.
%tfr=tfr1/(xrow/2);%the amplitude of tfr result has been pre-rectified.
Te=tfr.*IF;
%Te=Te/(xrow/2);
end
% %The following code is an alternative way to estimate IF.
% %In theroy, they are same.
% omega = zeros(round(N/2),tcol);
% for b=1:N
% omega(:,b) = (ft-1)'+real(va*1i*tfr2(ft,b)/2/pi./tfr1(ft,b));
% end
% for i=1:round(N/2)%frequency
% for j=1:N%time
%     if abs(tfr1(i,j))>0.8*E%default frequency resolution is 1Hz.
%         if abs(omega(i,j)-i)<0.5%default frequency resolution is 1Hz.
%         IF(i,j)=1;
%         end
%     end
% end
% end
