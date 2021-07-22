function [Data2]=SpectralMatching(Data,dt,TS,M,d,ss,maxgen,sizepop,n,o)
% Data    -- 输入时程
% dt      -- 输入时程的时间间隔
% TS      -- 匹配的周期矩阵（因JSDFYP函数周期为0.01的倍数，此处周期必须为0.01的倍数，JSDFYP_1无影响）
% M       -- 目标反应谱值（单列正值，输入前请与TS*对齐*!）
% d       -- 阻尼比（一般取0.05）
% ss      -- 迭代次数（一般取10）
% maxgen  -- 最大迭代次数（粒子群算法参数）（一般取5）
% sizepop -- 种群的规模，粒子个数（粒子群算法参数）（一般取12）
% n       -- 1:粒子群算法，2:伪穷举法，3:结合法（一般取1）
% o       -- 误差函数形式，决定适应度函数。1:总差值的均值，2:总比值的均值，3:单个差值的最大值，4:单个比值的最大值。（一般取1）

% 调用函数：Baseline_correction_JK_change,JSDFYP,JSDFYP_1,Cij,liziqun

%% 计算和设置参数
W=2*pi./TS;%周期矩阵对应的频率矩阵
TT=length(TS);
jgd=ss/5;%画图间隔
ss=ss+1;%循环次数加一
%% 画目标反应谱作为底图
% figure(1)
% set(gca,'xscale','log');%对数坐标系
% set(gca,'yscale','log');%对数坐标系
% hold on
% plot(TS,M,'-k','LineWidth',2)
% title('Sa-T')
% xlabel('T(s)')
% ylabel('Sa(cm/s2)')
% grid
%% 迭代计算，匹配反应谱
dT=dt;acc=Data;
[acc_correction,vel_correction,dis_correction]=Baseline_correction_JK_change(dT,acc);%基线校正
Data=acc_correction;
D(:,1)=Data;%D矩阵存储加速度时程
t=dt;GA=Data;
[PA,Ti,P]=JSDFYP_1(t,GA,TS,d);

% figure(1)
% hold on
% plot(TS,PA,'LineWidth',2)

Data=(max(M)/max(PA)).*Data;%加速度时程调幅
D(:,2)=Data;

for lp=1:ss
    t=dt;GA=Data;
    [PA,Ti,P]=JSDFYP_1(t,GA,TS,d);%计算加速度时程的反应谱、时间和极性
    m(:,lp)=PA;
    T(:,lp)=Ti;
    
    if lp==1||lp==2||(mod((lp-1),jgd)<0.1&&lp>2)
%         figure(1)
%         hold on
%         plot(TS,m(:,lp),'LineWidth',2)
    end
    
    for lp3=1:TT
        R1(lp3,1)=(M(lp3)-m(lp3,lp))*P(lp3);%R'=Q-R
    end
    
    m1(:,lp)=abs(m(:,lp)-M);%差值
    m2(:,lp)=m1(:,lp)./M;%比值
    m3(1,lp)=sum(m1(:,lp))./TT;%总差值的均值
    m3(2,lp)=sum(m2(:,lp))./TT;%总比值的均值
    m3(3,lp)=max(m1(:,lp));%单个差值的最大值
    m3(4,lp)=max(m2(:,lp));%单个比值的最大值
    
    disp(m3(o,lp))
    
    if lp==ss
        GA=Data;t=dt;t0=0.01;T0=10.0;
        [PA, BPA,T1]=JSDFYP(t,GA,t0,T0,d);
%         figure(2)
%         set(gca,'xscale','log');%对数坐标系
%         set(gca,'yscale','log');%对数坐标系
%         hold on
%         plot(TS,M,'-k',TS,m(:,1),'-r',T1,PA,'-g','LineWidth',2)
%         title('Sa-T')
%         xlabel('T(s)')
%         ylabel('Sa(cm/s2)')
%         legend('目标反应谱','原始地震动反应谱','合成地震动反应谱')
%         grid
        break
    end
    
    for lp4=1:TT
        for lp5=1:TT
            wi=W(lp4);wj=W(lp5);d=0.05;d1=sqrt(1-d^2);ti=T(lp4,lp);tj=T(lp5,lp);
            dtj=atan(d1/d)/(d1*wj);y=1.178*(wj/(2*pi))^(-0.93);
            C(lp4,lp5)=quadl(@Cij,0,ti,[],[],wi,wj,d1,d,ti,tj,dtj,y);%C
        end
    end
    
    B=C\R1;
    
    for lp6=1:length(Data)
        tt=lp6*dt;
        for lp7=1:TT
            wj=W(lp7);t=T(lp7,lp);dtj=atan(d1/d)/(d1*wj);y=1.178*(wj/(2*pi))^(-0.93);
            fj(lp6,lp7)=cos(d1.*wj.*(tt-t+dtj)).*exp(-((tt-t+dtj)./y).^2);
        end
    end
    
    if n==1
        [y1]=liziqun(maxgen,sizepop,o,M,TS,Data,fj,B,d,dt);%粒子群算法
    elseif n==2
        t1=clock;
        g=100;%穷举个数
        for i=1:g+1
            k=(i-1)/g;
            f(i)=fun(k,o,M,TS,Data,fj,B,d,dt);%伪穷举法
        end
        t2=clock;
        disp(etime(t2,t1))
        f1=f(1);ff1=1;
        for i=1:g+1
            if f(i)<f1
                f1=f(i);ff1=i;
            end
        end
        y1=(ff1-1)/g;
    elseif n==3
        %结合法
    end
    
    y2(lp)=y1;
    if y1==0%&&R2>1
        y1=0.001;
    end
    Data=Data+y1*fj*B;
    D(:,lp+2)=Data;
end
Data2=Data;

% figure(3)
% yyaxis left
% plot(1:ss-1,m3(1,2:end),'-r','LineWidth',2)
% hold on
% plot(1:ss-1,m3(3,2:end),'--r','LineWidth',2)
% set(gca,'yscale','log');%对数坐标系
% set(gca,'Ycolor','r')
% xlabel('迭代次数')
% ylabel('差值')
% axis([1 20 0.01 100])
% yyaxis right
% plot(1:ss-1,m3(2,2:end),'-k','LineWidth',2)
% plot(1:ss-1,m3(4,2:end),'--k','LineWidth',2)
% set(gca,'Ycolor','k')
% set(gca,'yscale','log');%对数坐标系
% ylabel('比值')
% legend('差值的总体均值','差值的个体最大值','比值的总体均值','比值的个体最大值')
% axis([1 20 0.0005 1])
% grid on
end
