function [acc_correction,vel_correction,dis_correction]=Baseline_correction_JK_change(dT,acc)
%%**********基线校正******************
%%&&&&&&&&&&采用一阶线性函数校正 a0+a1*T ********
T=dT*(1:length(acc))';
vel=cumsum(acc)*dT;
dis=cumsum(vel)*dT;
temp=cumsum(dis.*(3*T(end).*(T).^2-2*(T).^3))*dT;
a1=28/13*(1/(T(end)^2))*(2*vel(end)-15/(T(end).^5)*temp(end));
a0=vel(end)/T(end)-a1*T(end)/2;

acc_correction=acc-(a0+a1*T);
vel_correction=cumsum(acc_correction)*dT; %
dis_correction=cumsum(vel_correction)*dT;

end
