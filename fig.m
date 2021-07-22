clc;clear;
load('D(huige_M5_R20_Data60).mat')
figure(1)
set(gca,'xscale','log');%��������ϵ
set(gca,'yscale','log');%��������ϵ
hold on
plot(M(:,1),M(:,2),'-k',M(:,1),m(:,1),'-r',T1,PA,'-g','LineWidth',2)
title('Sa-T')

legend('Ŀ�귴Ӧ��','ԭʼ���𶯷�Ӧ��','�ϳɵ��𶯷�Ӧ��','Location','southwest')
axis([0.01 10 0.05 1000])
set(gca,'xscale','log');%��������ϵ
set(gca,'yscale','log');%��������ϵ
set(gca,'xTick',[0.01,0.03,0.1,0.3,1,3,10]);%����Ҫ��ʾ����̶�
set(gca,'yTick',[0.1 1 10 100 1000]);%����Ҫ��ʾ����̶�
grid on
box on
ylabel('\fontname{Times New Roman}Spectral Acceleration(gal)')
xlabel('\fontname{Times New Roman}Period(s)')
%text(0.55,35,'R_e_p_i=200km','color','k','FontSize',13)
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'FontSize',11.5)

Width=9;Height=8;%��λΪ���ף������������������ġ�����
ScreenSize=12; % ��Ļ��С����λΪӢ�磬��Ӧ��ע���ֵͨ��ָ�Խ��ߵĳ��ȣ�����ݹ��ɶ��������
ScreenSizeInCM=ScreenSize*2.45; %1Ӣ�����2.45���ף����Ȼ���
scrsz = get(0,'ScreenSize');  %�õ���Ļ�ֱ���
ScreenWidth=ScreenSizeInCM/sqrt(1+(scrsz(4)/scrsz(3))^2);%��Ļ����λΪ����
ScreenHeight=ScreenWidth*scrsz(4)/scrsz(3);%��Ļ�ߣ���λ����
WidthRatio=Width/ScreenWidth;%ͼ�ε������������Ļ��ȵı�ֵ
HeightRatio=Height/ScreenHeight;%ͼ�ε������߶�����Ļ�߶ȵı�ֵ
set(gcf,'Unit','Normalized','Position',[0.1 0.1 WidthRatio HeightRatio]);%���û�ͼ�Ĵ�
