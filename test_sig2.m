clc
clear
close all
%裁剪 32 99 76 191
fs=2^10;    %采样频率
dt=1/fs;    %时间精度
timestart=0;
timeend=1;
time=(0:(timeend-timestart)/dt-1)*dt+timestart;

Sig1 = exp(1j*(2*pi*(250*time))); %
Sig2 = exp(-0.5*time).*exp(1j*(500*pi*time + 100*sin(2*pi*time)));
Sig3 = 0.8*exp(0.5*time).*exp(1j*(500*pi*time - 100*sin(2*pi*time)));
sig = Sig1 + Sig2 + Sig3;
sig=sig(:);

N=length(sig);
fre = (fs/2)/(N/2):(fs/2)/(N/2):(fs/2);
df=(fs/2)/(N/2);




beta=3;
sigma=0.3;
hlength=128;%越大越抗干扰,越大对于我们的方法越有利
STFT = stft(sig, hlength,sigma);
[SST1,omega1,IF_est2] = SST(sig,hlength);
[SET1,IF_est3,~] = SET(sig,hlength);
[RM1,IF_est4] = RS(sig,hlength);
[LMSST1,omega5,IF_est5] = LMSST(sig,hlength,15);

%% ridge extraction and fitting
bw = fs/200;%% use the Fourier model to smooth the detected ridge curves；orderIF1 could be larger than the following orderIF
orderIF = 10;
num = 3; % the number of the components
delta = 10;
alpha = 1;
Nfrebin=length(fre);
window=hlength;
SampFreq=fs;
[fidexmult, tfdv] = extridge_mult(sig, fs, num, delta, orderIF,bw,Nfrebin,window,alpha,sigma);

%% ridge path regrouping
thrf = length(fre)/30;
[findex,interset] = RPRG(fidexmult,thrf);

%% signal decomposition
alpha = 0.5;
bw = SampFreq/40; % bandwidth of the ICCD
orderamp = round(bw*length(sig)/SampFreq);%% Fourier order for characterizing signal amplitudes
[extr_Sig,ampmatrix,IFfit] = ICCD(sig,SampFreq,fre(findex),orderIF,orderamp,alpha);% ICCD performs the joint-optimization scheme using all the obtained IFs

OSST = SQ(fre,STFT,IFfit);




figure()
imagesc(time,fre,abs(STFT));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',12,'FontName','Times New Roman');
set(gca,'FontSize',16,'FontName','Times New Roman')
set(gca,'linewidth',1);


figure()
imagesc(time,fre,abs(SST1));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',12,'FontName','Times New Roman');
set(gca,'FontSize',16,'FontName','Times New Roman')
set(gca,'linewidth',1);


figure()
imagesc(time,fre,abs(SET1));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',12,'FontName','Times New Roman');
set(gca,'FontSize',16,'FontName','Times New Roman')
set(gca,'linewidth',1);



figure()
imagesc(time,fre,abs(RM1));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',12,'FontName','Times New Roman');
set(gca,'FontSize',16,'FontName','Times New Roman')
set(gca,'linewidth',1);


figure()
imagesc(time,fre,abs(LMSST1));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',12,'FontName','Times New Roman');
set(gca,'FontSize',16,'FontName','Times New Roman')
set(gca,'linewidth',1);


figure()
imagesc(time,fre,abs(OSST));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',12,'FontName','Times New Roman');
set(gca,'FontSize',16,'FontName','Times New Roman')
set(gca,'linewidth',1);



Re1 =renyi(STFT);
Re2 =renyi(SST1);
Re3 =renyi(SET1);
Re4 = renyi(RM1);
Re5 = renyi(LMSST1);
Re6 = renyi(OSST);
RE=[Re1,Re2,Re3,Re4,Re5,Re6]

