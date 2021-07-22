function [PA,Ti,P]=JSDFYP_1(t,GA,TS,d)
% JSDFYP的改进版，计算周期点减少，增加计算了对应周期点最大加速度值的时间和正负

% t——地震动时间间隔
% GA——地震动时程
% TS——反应谱计算周期点
% d——反应谱阻尼比
%计算加速度反应谱
S=size(GA);
L=S(1)*S(2);
GA=reshape(GA,L,1);
TT=length(TS);
for lp=1:TT
    w=2*pi/TS(lp);
    a11=exp(-d*w*t)*(d*sin(w*sqrt(1-d^2)*t)/sqrt(1-d^2)+cos(w*sqrt(1-d^2)*t));
    a12=exp(-d*w*t)*sin(w*sqrt(1-d^2)*t)/(w*sqrt(1-d^2));
    a21=-w*exp(-d*w*t)*sin(w*sqrt(1-d^2)*t)/sqrt(1-d^2);
    a22=exp(-d*w*t)*(cos(w*sqrt(1-d^2)*t)-d*sin(w*sqrt(1-d^2)*t)/sqrt(1-d^2));
    b11=exp(-d*w*t)*(sin((w*(1-d^2)^(1/2))*t)/(w*sqrt(1-d^2))*((2*(d^2)-1)/(w^2*t)+d/w)+ cos(w*sqrt(1-d^2)*t)*(2*d/(w^3*t)+1/w^2))-2*d/(w^3*t);
    b12=-exp(-d*w*t)*((w*sqrt(1-d^2) )^(-1)*sin((w*(1-d^2)^(1/2))*t)*(2*(d^2)-1)/(w^2*t)+ cos(w*sqrt(1-d^2)*t)*2*d/(w^3*t))- 1/w^2+2*d/(w^3*t);
    b21=exp(-d*w*t)*((cos(w*sqrt(1-d^2)*t)-d*sqrt(1-d^2)^(-1)*sin((w*(1-d^2)^(1/2))*t))*((2*(d^2)-1)/(w^2*t)+d/w)-(w*sqrt(1-d^2)*sin((w*(1-d^2)^(1/2))*t)+d*w*cos(w*sqrt(1-d^2)*t))*(2*d/(w^3*t)+1/w^2))+1/(w^2*t);
    b22=-exp(-d*w*t)*((2*(d^2)-1)/(w^2*t)*(cos(w*sqrt(1-d^2)*t)-d*sqrt(1-d^2)^(-1)*sin((w*(1-d^2)^(1/2))*t))-2*d/(w^3*t)*(w*sqrt(1-d^2)*sin((w*(1-d^2)^(1/2))*t)+d*w*cos(w*sqrt(1-d^2)*t)))-1/(w^2*t);
    D(1)=0;
    V(1)=-GA(1)*t;
    A1(1)=2*d*w*GA(1)*t;
    for i=1:(L-1)
        D(i+1)=a11*D(i)+a12*V(i)+b11*GA(i)+b12*GA(i+1);
        V(i+1)=a21*D(i)+a22*V(i)+b21*GA(i)+b22*GA(i+1);
        A(i+1)=-GA(i+1)-2*d*w*V(i+1)-w^2*D(i+1);
        A1(i+1)=A(i+1)+GA(i+1);
    end
    
    [Y2,I2]=max(abs(A1));
    PA(lp)=Y2;
    Ti(lp)=I2*t;
    if A1(I2)>0
        P(lp)=1;
    else
        P(lp)=-1;
    end
end
end
