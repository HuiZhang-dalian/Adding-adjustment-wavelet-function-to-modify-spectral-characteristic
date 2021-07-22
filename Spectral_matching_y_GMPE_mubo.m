clc;clear
%% ���ò���
Mw=input('������𶯵��𼶣�');
R=input('����(Km)��');
% sta=input('���о���(������112.1��ʾ����112.1��)��');
% sto=input('����γ�ȣ�');
% epa=input('ָ���㾭�ȣ�');
% epo=input('ָ����γ�ȣ�');
% R=6371.004*acos(sind(sto)*sind(epo)+cosd(sto)*cosd(epo)*cosd(sta-epa));
TS=[0.01,0.02,0.03,0.04,0.05,0.08,0.1,0.12,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1,1.5,2,2.5,3,4,5,6,8,10];
W=2*pi./TS;
ss=20+1;TT=length(TS);jgd=2;
maxgen=5;%����������
sizepop=12;%��Ⱥ�Ĺ�ģ��100������
%% ����Ŀ���ײ���ͼ
[M]=target(Mw,R);

figure(1)
set(gca,'xscale','log');%��������ϵ
set(gca,'yscale','log');%��������ϵ
hold on
plot(M(:,1),M(:,2),'-k','LineWidth',2)
title('Sa-T')
xlabel('T(s)')
ylabel('Sa(cm/s2)')
grid
%% �������㣬ƥ�䷴Ӧ��
load('Data2.mat');%�����ʼ��������
D(:,1)=Data2;
%load('DATA.mat');%�����ʼ��������
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
        set(gca,'xscale','log');%��������ϵ
        set(gca,'yscale','log');%��������ϵ
        hold on
        plot(M(:,1),m(:,lp),'LineWidth',2)
        title('Sa-T')
        xlabel('T(s)')
        ylabel('Sa(cm/s2)')
        grid
    elseif mod(lp,jgd)<0.1
        figure(1)
        set(gca,'xscale','log');%��������ϵ
        set(gca,'yscale','log');%��������ϵ
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
        set(gca,'xscale','log');%��������ϵ
        set(gca,'yscale','log');%��������ϵ
        hold on
        plot(M(:,1),M(:,2),'-k',M(:,1),m(:,1),'-r',T1,PA,'-g','LineWidth',2)
        title('Sa-T')
        xlabel('T(s)')
        ylabel('Sa(cm/s2)')
        legend('Ŀ�귴Ӧ��','ԭʼ���𶯷�Ӧ��','�ϳɵ��𶯷�Ӧ��')
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
        
        %% ����Ⱥ�㷨�����Ž�
        c1=2.05;%ÿ�����ӵĸ���ѧϰ���ӣ����ٶȳ���
        c2=2.05;%ÿ�����ӵ����ѧϰ���ӣ����ٶȳ���
        Vmax=3;%���ӵ��������ٶ�
        Vmin=-3;%���ӵ��������ٶ�
        popmax=1; %����xȡֵ-10-10
        popmin=0;
        ws=0.9; %wȨֵ���ٶ�����
        we=0.4; %wȨֵ���ٶ�����
        for i=1:sizepop-2
            % �������һ����Ⱥ
            pop(i,:)=rand(1,1); %��ʼ��Ⱥ�����ʵ�ֶ��ڷ�Χ-100-100��rand����Ϊ0��1֮��
            V(i,:)=3*rands(1,1); %��ʼ�ٶȣ��붨����ٶ����ֵ����Сֵһ�� -3-3,
            %������Ӧ�Ⱥ���
            fitness(i)=fun(pop(i,1));
        end
        pop(sizepop-1,:)=0; V(sizepop-1,:)=3*rands(1,1); fitness(sizepop-1)=fun(pop(sizepop-1,1));
        pop(sizepop,:)=1; V(sizepop,:)=3*rands(1,1); fitness(sizepop)=fun(pop(sizepop,1));
        %% ���弫ֵ��Ⱥ�弫ֵ
        [bestfiness,bestindex]=min(fitness); %beatfiness-��ֵ��beatindex-λ��
        zbest=pop(bestindex,:);%ȫ������
        gbest=pop;%��������
        fitnessgbest=fitness; %���������Ӧ��ֵ����ȫ�����Ž�
        fitnesszbest=bestfiness; %ȫ�������Ӧ��ֵ
        
        %% ����Ѱ��
        for i=1:maxgen
            w=ws-(ws-we)*(i/maxgen);
            %����50��
            for j=1:sizepop
                %�ٶȸ���
                V(j,:)=V(j,:)+c1*rand*(gbest(j,:)-pop(j,:))+c2*rand*(zbest-pop(j,:)); %�ٶȸ��¹�ʽ
                %V(j,:)=w*V(j,:)+c1*rand*(gbest(j,:)-pop(j,:))+c2*rand*(zbest-pop(j,:)); %�ٶȸ��¹�ʽ
                V(j,find(V(j,:)>Vmax))=Vmax; %�߽�Լ��
                V(j,find(V(j,:)<Vmin))=Vmin;%�߽�Լ��
                %��Ⱥ����
                pop(j,:)=pop(j,:)+V(j,:);
                pop(j,find(pop(j,:)>popmax))=popmax;%�߽�Լ��
                pop(j,find(pop(j,:)<popmin))=popmin;%�߽�Լ��
                %�ʶȸ���
                fitness(j)=fun(pop(j,1));
            end
            
            for j=1:sizepop
                %��������
                if fitness(j)<fitnessgbest(j)
                    gbest(j,:)=pop(j,:);
                    fitnessgbest(j)=fitness(j);
                    y1=gbest(j,:);
                end
                %Ⱥ������
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
    m1(:,i)=abs(m(:,i)-M(:,2));%��ֵ
    m2(:,i)=m1(:,i)./M(:,2);%��ֵ
end

m3(1,:)=sum(m1)./30;%��ֵ�������ֵ
m3(2,:)=sum(m2)./30;%��ֵ�������ֵ
m3(3,:)=max(m1);%��ֵ�ĸ������ֵ
m3(4,:)=max(m2);%��ֵ�ĸ������ֵ

figure(3)
yyaxis left
plot(1:ss-1,m3(1,2:end),'-r','LineWidth',2)
hold on
plot(1:ss-1,m3(3,2:end),'--r','LineWidth',2)
set(gca,'yscale','log');%��������ϵ
set(gca,'Ycolor','r')
xlabel('��������')
ylabel('��ֵ')
axis([1 20 0.01 100])
yyaxis right
plot(1:ss-1,m3(2,2:end),'-k','LineWidth',2)
plot(1:ss-1,m3(4,2:end),'--k','LineWidth',2)
set(gca,'Ycolor','k')
set(gca,'yscale','log');%��������ϵ
ylabel('��ֵ')
legend('��ֵ�������ֵ','��ֵ�ĸ������ֵ','��ֵ�������ֵ','��ֵ�ĸ������ֵ')
axis([1 20 0.0005 1])
grid on

