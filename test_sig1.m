clear;
close all;
clc;
%�ü� 35 99 78 191
%�ü� 15 89 33 176
fs=2^7;    %����Ƶ��
dt=1/fs;    %ʱ�侫��
timestart=0;
timeend=1;
time=(0:(timeend-timestart)/dt-1)*dt+timestart;
Sig1=sin(2*pi*(40*time + 1*sin(4*pi*time)));
Sig2 = sin(2*pi*(10*time + 10*(time-0.5).^3));
sig=Sig1+Sig2;
sig=sig(:);
N=length(sig);
fre = (fs/2)/(N/2):(fs/2)/(N/2):(fs/2);
df=(fs/2)/(N/2);
alpha=(fs/2)/(N/2);


beta=3;
sigma=0.2;
hlength=50;%Խ��Խ������,Խ��������ǵķ���Խ����
STFT = stft(sig, hlength,sigma);
[SST1,omega1,IF_est2] = SST(sig,hlength);
[SET1,IF_est3,~] = SET(sig,hlength);
[RM1,IF_est4] = RS(sig,hlength);
[LMSST1,omega5,IF_est5] = LMSST(sig,hlength,15);
[OSST,OSST2,Dtfr2,omega6,IF_est6,T2] = SQ2(time,STFT,fre,0.6);
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
Re6 = renyi(OSST2);
RE=[Re1,Re2,Re3,Re4,Re5,Re6]


width=800;%��ȣ�������
height=fix(600*0.8);%�߶�
left=200;%����Ļ���½�ˮƽ����
bottem=100;%����Ļ���½Ǵ�ֱ����

figure()
set(gcf,'position',[left,bottem,width,height])  
subplot(231)
imagesc(time,fre,abs(STFT));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',10,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',10,'FontName','Times New Roman');
set(gca,'FontSize',12,'FontName','Times New Roman')
set(gca,'linewidth',1);
title('STFT')
subplot(232)
imagesc(time,fre,abs(SST1));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',10,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',10,'FontName','Times New Roman');
set(gca,'FontSize',12,'FontName','Times New Roman')
set(gca,'linewidth',1);
title('SST')

subplot(233)
imagesc(time,fre,abs(SET1));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',10,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',10,'FontName','Times New Roman');
set(gca,'FontSize',12,'FontName','Times New Roman')
set(gca,'linewidth',1);
title('SET')


subplot(234)
imagesc(time,fre,abs(RM1));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',10,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',10,'FontName','Times New Roman');
set(gca,'FontSize',12,'FontName','Times New Roman')
set(gca,'linewidth',1);
title('RM')

subplot(235)
imagesc(time,fre,abs(LMSST1));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',10,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',10,'FontName','Times New Roman');
set(gca,'FontSize',12,'FontName','Times New Roman')
set(gca,'linewidth',1);
title('LMSST')

subplot(236)
imagesc(time,fre,abs(OSST));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',10,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',10,'FontName','Times New Roman');
set(gca,'FontSize',12,'FontName','Times New Roman')
set(gca,'linewidth',1);
title('Ours')
