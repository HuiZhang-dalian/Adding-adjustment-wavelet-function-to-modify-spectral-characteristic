clc;clear;
load('D(huige_M5_R20_Data60).mat')
figure(1)
set(gca,'xscale','log');%对数坐标系
set(gca,'yscale','log');%对数坐标系
hold on
plot(M(:,1),M(:,2),'-k',M(:,1),m(:,1),'-r',T1,PA,'-g','LineWidth',2)
title('Sa-T')

legend('目标反应谱','原始地震动反应谱','合成地震动反应谱','Location','southwest')
axis([0.01 10 0.05 1000])
set(gca,'xscale','log');%对数坐标系
set(gca,'yscale','log');%对数坐标系
set(gca,'xTick',[0.01,0.03,0.1,0.3,1,3,10]);%设置要显示坐标刻度
set(gca,'yTick',[0.1 1 10 100 1000]);%设置要显示坐标刻度
grid on
box on
ylabel('\fontname{Times New Roman}Spectral Acceleration(gal)')
xlabel('\fontname{Times New Roman}Period(s)')
%text(0.55,35,'R_e_p_i=200km','color','k','FontSize',13)
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'FontSize',11.5)

Width=9;Height=8;%单位为厘米！！！这里根据需求更改。。。
ScreenSize=12; % 屏幕大小，单位为英寸，且应该注意该值通常指对角线的长度，需根据勾股定理计算宽高
ScreenSizeInCM=ScreenSize*2.45; %1英寸等于2.45厘米，长度换算
scrsz = get(0,'ScreenSize');  %得到屏幕分辨率
ScreenWidth=ScreenSizeInCM/sqrt(1+(scrsz(4)/scrsz(3))^2);%屏幕宽，单位为厘米
ScreenHeight=ScreenWidth*scrsz(4)/scrsz(3);%屏幕高，单位厘米
WidthRatio=Width/ScreenWidth;%图形的期望宽度与屏幕宽度的比值
HeightRatio=Height/ScreenHeight;%图形的期望高度与屏幕高度的比值
set(gcf,'Unit','Normalized','Position',[0.1 0.1 WidthRatio HeightRatio]);%设置绘图的大
