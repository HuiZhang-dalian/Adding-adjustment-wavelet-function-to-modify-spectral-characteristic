%clear;clc;
t=0.01;GA=Data;t0=0.01;T0=10;d=0.05;
[PA, BPA,T1]=JSDFYP(t,GA,t0,T0,d);

figure(1)
plot(T1,PA,'LineWidth',2)
set(gca,'xscale','log');%对数坐标系
set(gca,'yscale','log');%对数坐标系
grid