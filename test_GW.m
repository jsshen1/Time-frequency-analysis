clear;
clc;
close all;
%²Ã¼ô 56 134 117 267
load('GWdata.mat')
load('GWdata_relativity.mat')

t=0:GWdata(2,1)-GWdata(1,1):GWdata(end,1)-GWdata(1,1);
f=1/(GWdata(2,1)-GWdata(1,1));
n=length(GWdata(:,2));
time=GWdata(:,1);
fre=(f/2)/(n/2):(f/2)/(n/2):(f/2);
sig=GWdata(:,2);

sigma=0.32;
hlength=600;%600
tic
STFT= stft(sig,hlength,sigma);
[SST1,~] = SST(sig,hlength);
[SET1,~] = SET(sig,hlength);
[RM1,~] = RS(sig,hlength);
[LMSST1,~] = LMSST(sig,hlength,23);
[OSST,OSST0,~,~,~,~] = SQ2(time,STFT,fre,0.1);
toc



tfr2=STFT;

L=80;
fre2=fre(1:L);
STFT2=STFT(1:L,:);
SST2=SST1(1:L,:);
SET2=SET1(1:L,:);
RM2=RM1(1:L,:);
LMSST2=LMSST1(1:L,:);
OSST2=OSST(1:L,:);


width=400;%¿í¶È£¬ÏñËØÊý
height=fix(400*0.3);%¸ß¶È
left=200;%¾àÆÁÄ»×óÏÂ½ÇË®Æ½¾àÀë
bottem=100;%¾àÆÁÄ»×óÏÂ½Ç´¹Ö±¾àÀë

figure()
set(gcf,'position',[left,bottem,width,height])  
imagesc(time,fre2,abs(STFT2));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',8,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',8,'FontName','Times New Roman');
set(gca,'FontSize',8,'FontName','Times New Roman')
set(gca,'linewidth',1);

figure()
set(gcf,'position',[left,bottem,width,height])  
imagesc(time,fre2,abs(SST2));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',8,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',8,'FontName','Times New Roman');
set(gca,'FontSize',8,'FontName','Times New Roman')
set(gca,'linewidth',1);

figure()
set(gcf,'position',[left,bottem,width,height])  
imagesc(time,fre2,abs(SET2));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',8,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',8,'FontName','Times New Roman');
set(gca,'FontSize',8,'FontName','Times New Roman')
set(gca,'linewidth',1);


figure() 
set(gcf,'position',[left,bottem,width,height])  
imagesc(time,fre2,abs(RM2));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',8,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',8,'FontName','Times New Roman');
set(gca,'FontSize',8,'FontName','Times New Roman')
set(gca,'linewidth',1);

figure()
set(gcf,'position',[left,bottem,width,height])  
imagesc(time,fre2,abs(LMSST2));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',8,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',8,'FontName','Times New Roman');
set(gca,'FontSize',8,'FontName','Times New Roman')
set(gca,'linewidth',1);

figure()
set(gcf,'position',[left,bottem,width,height])  
imagesc(time,fre2,abs(OSST2));
set(gca,'ydir','normal');
xlabel('Time(s)','FontSize',8,'FontName','Times New Roman');
ylabel('Frequency(Hz)','FontSize',8,'FontName','Times New Roman');
set(gca,'FontSize',8,'FontName','Times New Roman')
set(gca,'linewidth',1);


Re1 =renyi(STFT);
Re2 =renyi(SST1);
Re3 =renyi(SET1);
Re4 = renyi(RM1);
Re5 = renyi(LMSST1);
Re6 = renyi(OSST0);

RE=[Re1,Re2,Re3,Re4,Re5,Re6]


