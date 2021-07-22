clc;clear
%% 设置参数
Mw=input('输入地震动的震级：');
R=input('距离(Km)：');
% sta=input('震中经度(如输入112.1表示东经112.1度)：');
% sto=input('震中纬度：');
% epa=input('指定点经度：');
% epo=input('指定点纬度：');
% R=6371.004*acos(sind(sto)*sind(epo)+cosd(sto)*cosd(epo)*cosd(sta-epa));
TS=[0.01,0.02,0.03,0.04,0.05,0.08,0.1,0.12,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1,1.5,2,2.5,3,4,5,6,8,10];
W=2*pi./TS;
ss=20+1;TT=length(TS);jgd=2;
maxgen=5;%最大迭代次数
sizepop=12;%种群的规模，100个粒子
%% 计算目标谱并画图
[M]=target(Mw,R);

figure(1)
set(gca,'xscale','log');%对数坐标系
set(gca,'yscale','log');%对数坐标系
hold on
plot(M(:,1),M(:,2),'-k','LineWidth',2)
title('Sa-T')
xlabel('T(s)')
ylabel('Sa(cm/s2)')
grid
%% 迭代计算，匹配反应谱
load('Data2.mat');%导入初始地震动数据
D(:,1)=Data2;
%load('DATA.mat');%导入初始地震动数据
%Data2=Data(:,40);
dt2=0.005;
dT=dt2;acc=Data2;
[ acc_correction,vel_correction,dis_correction] = Baseline_correction_JK_change( dT,acc);
Data2=acc_correction;

b=0.05;
for lp=1:ss
    k2=1;
    for j2=1:1:TT
        T2=TS(j2);
        w=2*pi/T2;d1(1)=0;v1(1)=0;a1(1)=0;
        for j1=2:size(Data2,1)+1
            d1(j1)=(-Data2(j1-1)+2*b*w*(3*d1(j1-1)/dt2+2*v1(j1-1)+dt2*a1(j1-1)/2)+...
                (6*d1(j1-1)/(dt2*dt2)+6*v1(j1-1)/dt2+2*a1(j1-1)))/(w^2+6*b*w/dt2+6/(dt2*dt2));
            v1(j1)=3*(d1(j1)-d1(j1-1))/dt2-2*v1(j1-1)-dt2*a1(j1-1)/2;
            a1(j1)=-Data2(j1-1)-2*b*w*v1(j1)-w*w*d1(j1);
            A1(j1)=a1(j1)+Data2(j1-1);
        end
        [Y2,I2]=max(abs(A1));
        T(k2,lp)=(I2-1)*dt2;
        m(k2,lp)=Y2;
        if A1(I2)>0
            P(k2)=1;
        else
            P(k2)=-1;
        end
        k2=k2+1;
    end
    
    if lp==1||lp==2
        figure(1)
        set(gca,'xscale','log');%对数坐标系
        set(gca,'yscale','log');%对数坐标系
        hold on
        plot(M(:,1),m(:,lp),'LineWidth',2)
        title('Sa-T')
        xlabel('T(s)')
        ylabel('Sa(cm/s2)')
        grid
    elseif mod(lp,jgd)<0.1
        figure(1)
        set(gca,'xscale','log');%对数坐标系
        set(gca,'yscale','log');%对数坐标系
        hold on
        plot(M(:,1),m(:,lp),'LineWidth',2)
        title('Sa-T')
        xlabel('T(s)')
        ylabel('Sa(cm/s2)')
        grid
    end
    
    if lp==ss
        GA=Data2;t=0.005;t0=0.01;T0=10.0;d=0.05;
        [PA, BPA,T1]=JSDFYP(t,GA,t0,T0,d);
        %save('D(huige_M4_R50_Data40).mat','M','m','T1','PA','Data')
        figure(2)
        set(gca,'xscale','log');%对数坐标系
        set(gca,'yscale','log');%对数坐标系
        hold on
        plot(M(:,1),M(:,2),'-k',M(:,1),m(:,1),'-r',T1,PA,'-g','LineWidth',2)
        title('Sa-T')
        xlabel('T(s)')
        ylabel('Sa(cm/s2)')
        legend('目标反应谱','原始地震动反应谱','合成地震动反应谱')
        grid
        break
    end
    
    if lp==1
        Data2=(max(abs(M(:,2)))/max(abs(m(:,lp)))).*Data2;
        D(:,2)=Data2;
    else
        R2=0;
        for lp3=1:TT
            R1(lp3,1)=(M(lp3,2)-m(lp3,lp))*P(lp3);
            R2=R2+abs(R1(lp3,1));
        end%R'=Q-R
        
        disp(R2)
        %R1(lp,1)=max(abs(R));
        %disp(max(abs(R)))
        
        for lp4=1:1:TT
            for lp5=1:1:TT
                wi=W(lp4);wj=W(lp5);b=0.05;b1=sqrt(1-b^2);ti=T(lp4,lp);tj=T(lp5,lp);
                dtj=atan(b1/b)/(b1*wj);y=1.178*(wj/(2*pi))^(-0.93);
                C(lp4,lp5)=quadl(@Cij,0,ti,[],[],wi,wj,b1,b,ti,tj,dtj,y);
            end
        end%Cij
        
        B=C\R1;
        
        for lp6=1:1:length(Data2)
            tt=lp6*dt2;
            for lp7=1:1:TT
                wj=W(lp7);t=T(lp7,lp);dtj=atan(b1/b)/(b1*wj);y=1.178*(wj/(2*pi))^(-0.93);
                fj(lp6,lp7)=cos(b1.*wj.*(tt-t+dtj)).*exp(-((tt-t+dtj)./y).^2);
            end
        end
        
        %% 粒子群算法求最优解
        c1=2.05;%每个粒子的个体学习因子，加速度常数
        c2=2.05;%每个粒子的社会学习因子，加速度常数
        Vmax=3;%粒子的最大飞翔速度
        Vmin=-3;%粒子的最大飞翔速度
        popmax=1; %根据x取值-10-10
        popmin=0;
        ws=0.9; %w权值，速度相结合
        we=0.4; %w权值，速度相结合
        for i=1:sizepop-2
            % 随机产生一个种群
            pop(i,:)=rand(1,1); %初始种群，如何实现定于范围-100-100，rand函数为0到1之间
            V(i,:)=3*rands(1,1); %初始速度，与定义的速度最大值和最小值一致 -3-3,
            %计算适应度函数
            fitness(i)=fun(pop(i,1));
        end
        pop(sizepop-1,:)=0; V(sizepop-1,:)=3*rands(1,1); fitness(sizepop-1)=fun(pop(sizepop-1,1));
        pop(sizepop,:)=1; V(sizepop,:)=3*rands(1,1); fitness(sizepop)=fun(pop(sizepop,1));
        %% 个体极值和群体极值
        [bestfiness,bestindex]=min(fitness); %beatfiness-数值，beatindex-位置
        zbest=pop(bestindex,:);%全局最优
        gbest=pop;%个体最优
        fitnessgbest=fitness; %个体最佳适应度值，即全局最优解
        fitnesszbest=bestfiness; %全局最佳适应度值
        
        %% 迭代寻优
        for i=1:maxgen
            w=ws-(ws-we)*(i/maxgen);
            %迭代50次
            for j=1:sizepop
                %速度更新
                V(j,:)=V(j,:)+c1*rand*(gbest(j,:)-pop(j,:))+c2*rand*(zbest-pop(j,:)); %速度更新公式
                %V(j,:)=w*V(j,:)+c1*rand*(gbest(j,:)-pop(j,:))+c2*rand*(zbest-pop(j,:)); %速度更新公式
                V(j,find(V(j,:)>Vmax))=Vmax; %边界约束
                V(j,find(V(j,:)<Vmin))=Vmin;%边界约束
                %种群更新
                pop(j,:)=pop(j,:)+V(j,:);
                pop(j,find(pop(j,:)>popmax))=popmax;%边界约束
                pop(j,find(pop(j,:)<popmin))=popmin;%边界约束
                %适度更新
                fitness(j)=fun(pop(j,1));
            end
            
            for j=1:sizepop
                %个体最优
                if fitness(j)<fitnessgbest(j)
                    gbest(j,:)=pop(j,:);
                    fitnessgbest(j)=fitness(j);
                    y1=gbest(j,:);
                end
                %群体最优
                if fitness(j)<fitnesszbest
                    zbest=pop(j,:);
                    fitnesszbest=fitness(j);
                    y1=zbest;
                end
            end
        end
        y2(lp)=y1;
          if y1==0&&R2>1
              y1=0.001;
          end
        Data2=Data2+y1*fj*B;
        D(:,lp+1)=Data2;
    end
end
Data=Data2;
save('Data1.mat','Data')

for i=2:ss
    m1(:,i)=abs(m(:,i)-M(:,2));%差值
    m2(:,i)=m1(:,i)./M(:,2);%比值
end

m3(1,:)=sum(m1)./30;%差值的总体均值
m3(2,:)=sum(m2)./30;%比值的总体均值
m3(3,:)=max(m1);%差值的个体最大值
m3(4,:)=max(m2);%比值的个体最大值

figure(3)
yyaxis left
plot(1:ss-1,m3(1,2:end),'-r','LineWidth',2)
hold on
plot(1:ss-1,m3(3,2:end),'--r','LineWidth',2)
set(gca,'yscale','log');%对数坐标系
set(gca,'Ycolor','r')
xlabel('迭代次数')
ylabel('差值')
axis([1 20 0.01 100])
yyaxis right
plot(1:ss-1,m3(2,2:end),'-k','LineWidth',2)
plot(1:ss-1,m3(4,2:end),'--k','LineWidth',2)
set(gca,'Ycolor','k')
set(gca,'yscale','log');%对数坐标系
ylabel('比值')
legend('差值的总体均值','差值的个体最大值','比值的总体均值','比值的个体最大值')
axis([1 20 0.0005 1])
grid on

