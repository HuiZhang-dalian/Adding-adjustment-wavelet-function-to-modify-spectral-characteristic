clear;clc;
Mw=input('输入地震动的震级：');
R=input('距离(Km)：');
% sta=input('震中经度(如输入112.1表示东经112.1度)：');
% sto=input('震中纬度：');
% epa=input('指定点经度：');
% epo=input('指定点纬度：');
% R=6371.004*acos(sind(sto)*sind(epo)+cosd(sto)*cosd(epo)*cosd(sta-epa));
load D.mat
Data=D(:,1);
dt=0.005;
TS=[0.01,0.02,0.03,0.04,0.05,0.08,0.1,0.12,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1,1.5,2,2.5,3,4,5,6,8,10];
[M]=target(Mw,R);M=M(:,2);
d=0.05;ss=10;maxgen=5;sizepop=12;n=2;o=1;
[Data2]=SpectralMatching(Data,dt,TS,M,d,ss,maxgen,sizepop,n,o);
