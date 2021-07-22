function f=fun(k,o,M,TS,Data,fj,B,d,dt)
TT=length(TS);
t=dt;%反应谱的时间间隔
GA=Data+k.*fj*B;
t0=0.01;T0=10.0;%反应谱的参数
PA=JSDFYP(t,GA,t0,T0,d);
for i=1:TT
    PA1(i,1)=PA(100*TS(i));
end

M1(:,1)=abs(PA1-M);%差值
M2(:,1)=M1(:,1)./M;%比值
M3(1,1)=sum(M1(:,1))./TT;%总差值的均值
M3(2,1)=sum(M2(:,1))./TT;%总比值的均值
M3(3,1)=max(M1(:,1));%单个差值的最大值
M3(4,1)=max(M2(:,1));%单个比值的最大值
f=M3(o,1);
end