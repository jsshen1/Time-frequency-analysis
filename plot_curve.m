clear;
clc;
close all;

f=0:0.001:10;
f=f';
miu=5;
g=1/sqrt(2*pi)*exp(-0.5*(f-miu).*(f-miu));
S=[f g];
figure()
plot(g);
set(gca,'ydir','normal');
% xlabel('Frequency(Hz)','FontSize',12,'FontName','Times New Roman');
% ylabel('Amp','FontSize',12,'FontName','Times New Roman');
set(gca,'FontSize',16,'FontName','Times New Roman')
set(gca,'linewidth',1);
title('STFT')