function [y1]=liziqun(maxgen,sizepop,o,M,TS,Data,fj,B,d,dt)
% 调用函数：fun
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
    fitness(i)=fun(pop(i,1),o,M,TS,Data,fj,B,d,dt);
end
pop(sizepop-1,:)=0; V(sizepop-1,:)=3*rands(1,1); fitness(sizepop-1)=fun(pop(sizepop-1,1),o,M,TS,Data,fj,B,d,dt);
pop(sizepop,:)=1; V(sizepop,:)=3*rands(1,1); fitness(sizepop)=fun(pop(sizepop,1),o,M,TS,Data,fj,B,d,dt);
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
        fitness(j)=fun(pop(j,1),o,M,TS,Data,fj,B,d,dt);
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
end
