function [M]=target(Mw,R)
load('COEF.mat');
TS=evalin('base', 'TS');
M(:,1)=TS;
for k=1:length(TS)
    M(k,2)=exp(COEF(k+1,1)+COEF(k+1,2).*(Mw-6)+COEF(k+1,3).*(Mw-6).^2+COEF(k+1,4).*log(sqrt(R.^2+COEF(k+1,5).^2)));
end