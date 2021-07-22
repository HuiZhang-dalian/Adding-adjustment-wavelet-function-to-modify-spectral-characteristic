clc;clear
T=[0.0100000000000000,0.0200000000000000,0.0220000000000000,0.0250000000000000,0.0290000000000000,0.0300000000000000,0.0320000000000000,0.0350000000000000,0.0360000000000000,0.0400000000000000,0.0420000000000000,0.0440000000000000,0.0450000000000000,0.0460000000000000,0.0480000000000000,0.0500000000000000,0.0550000000000000,0.0600000000000000,0.0650000000000000,0.0670000000000000,0.0700000000000000,0.0750000000000000,0.0800000000000000,0.0850000000000000,0.0900000000000000,0.0950000000000000,0.100000000000000,0.110000000000000,0.120000000000000,0.130000000000000,0.133000000000000,0.140000000000000,0.150000000000000,0.160000000000000,0.170000000000000,0.180000000000000,0.190000000000000,0.200000000000000,0.220000000000000,0.240000000000000,0.250000000000000,0.260000000000000,0.280000000000000,0.290000000000000,0.300000000000000,0.320000000000000,0.340000000000000,0.350000000000000,0.360000000000000,0.380000000000000,0.400000000000000,0.420000000000000,0.440000000000000,0.450000000000000,0.460000000000000,0.480000000000000,0.500000000000000,0.550000000000000,0.600000000000000,0.650000000000000,0.667000000000000,0.700000000000000,0.750000000000000,0.800000000000000,0.850000000000000,0.900000000000000,0.950000000000000,1,1.10000000000000,1.20000000000000,1.30000000000000,1.40000000000000,1.50000000000000,1.60000000000000,1.70000000000000,1.80000000000000,1.90000000000000,2,2.20000000000000,2.40000000000000,2.50000000000000,2.60000000000000,2.80000000000000,3,3.20000000000000,3.40000000000000,3.50000000000000,3.60000000000000,3.80000000000000,4,4.20000000000000,4.40000000000000,4.60000000000000,4.80000000000000,5,5.50000000000000,6,6.50000000000000,7,7.50000000000000,8,8.50000000000000,9,9.50000000000000,10]
Mw=6;
Rjb=10
Fault_Type=0
region=1
z1=999
Vs30=360
T1=1 %Tcond
e=1 %
k=1 % one standard
[median, sigma, period1] = gmpe_bssa_2014(Mw, T, Rjb, Fault_Type, region, z1, Vs30);
nn=length(T);
for i = 1:nn
    rho(i) = baker_jayaram_correlation(T(i), T1); %correlation relationship 如果采用RJK16（最大周期不超2.0s）
end

for i=1:nn
meanReq(i,1)=T(i);
meanReq(i,2) =exp(log(median(i))+sigma(i).*e.*rho(i)); % 单位（gal）
end

for i=1:nn
meanReq(i,3)=exp(log(meanReq(i,2))+k*sigma(i)*sqrt(1-rho(i)^2))
meanReq(i,4)=exp(log(meanReq(i,2))-k*sigma(i)*sqrt(1-rho(i)^2))
end

result=meanReq 
figure 
hold on
Loca=find(T1==T)
plot(result(Loca,1), result(Loca,2),'o','linewidth',1.5,'color','r','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','b');

h0=plot(result(:,1),result(:,2),'linewidth',2,'color','b','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k');



h1=plot(result(:,1),result(:,3),'linewidth',2,'color','r','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k');
plot(result(:,1),result(:,4),'linewidth',2,'color','r','MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor','b');

% h2=plot(selectionParams.TgtPer, exp(targetSa.meanReq + 1.96*sqrt(diag(targetSa.covReq))'), '--r','linewidth',2,'color','r','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k');
% plot(selectionParams.TgtPer, exp(targetSa.meanReq - 1.96*sqrt(diag(targetSa.covReq))'), '--r','linewidth',2,'color','r','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k');
% h3=plot(knownPer,SaKnown(IMs.recID(1),:).*repmat(IMs.scaleFac(1),1,size(SaKnown,2)),'linewidth',1.0,'color','k','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k');
% plot(knownPer,SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2)),'linewidth',1.0,'color','k','MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor','k');
axis([0.01 10 1e-2 10])
xlabel('\fontname{Times New Roman}\itT\rm (s)')
ylabel('\fontname{Times New Roman}S_a (g)')
text(0.5,0.1,'\fontname{Times New Roman}\itT\rm_c_o_n_d')

figureFontSize=13
box on
legend([h0,h1],'\fontname{Times New Roman}Median','\fontname{Times New Roman}2.5 and 97.5 percentile')
set(findall(gcf,'-property','FontSize'),'FontSize', figureFontSize)
set(gca,'xscale','log');%对数坐标系
set(gca,'yscale','log');%对数坐标系
set(gca,'LineWidth',1.2)



Width=6.5;Height=5.1;%单位为厘米！！！这里根据需求更改。。。。。
ScreenSize=10.6; % 屏幕大小，单位为英寸，且应该注意该值通常指对角线的长度，需根据勾股定理计算宽高
ScreenSizeInCM=ScreenSize*2.45; %1英寸等于2.45厘米，长度换算
scrsz = get(0,'ScreenSize');  %得到屏幕分辨率
ScreenWidth=ScreenSizeInCM/sqrt(1+(scrsz(4)/scrsz(3))^2);%屏幕宽，单位为厘米
ScreenHeight=ScreenWidth*scrsz(4)/scrsz(3);%屏幕高，单位厘米
WidthRatio=Width/ScreenWidth;%图形的期望宽度与屏幕宽度的比值
HeightRatio=Height/ScreenHeight;%图形的期望高度与屏幕高度的比值
set(gcf,'Unit','Normalized','Position',[0.1 0.1 WidthRatio HeightRatio]);%设置绘图的大

