function f=fun(k,o,M,TS,Data,fj,B,d,dt)
TT=length(TS);
t=dt;%��Ӧ�׵�ʱ����
GA=Data+k.*fj*B;
t0=0.01;T0=10.0;%��Ӧ�׵Ĳ���
PA=JSDFYP(t,GA,t0,T0,d);
for i=1:TT
    PA1(i,1)=PA(100*TS(i));
end

M1(:,1)=abs(PA1-M);%��ֵ
M2(:,1)=M1(:,1)./M;%��ֵ
M3(1,1)=sum(M1(:,1))./TT;%�ܲ�ֵ�ľ�ֵ
M3(2,1)=sum(M2(:,1))./TT;%�ܱ�ֵ�ľ�ֵ
M3(3,1)=max(M1(:,1));%������ֵ�����ֵ
M3(4,1)=max(M2(:,1));%������ֵ�����ֵ
f=M3(o,1);
end