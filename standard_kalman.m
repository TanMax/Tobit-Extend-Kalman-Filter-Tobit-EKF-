%以温度测量为例 2020/2/6 谭承旦
clc;clear all,clear figure;
N=3000;%采样点的个数
CON=25;%室内温度理论值
Xexpect=CON*ones(1,N);%期望温度
X=zeros(1,N);%初始化真实温度值  
Xkf=zeros(1,N);%kalman滤波器估计值
Z=zeros(1,N);%温度计测量值
P=zeros(1,N);%初值化协方差 

%赋初值
X(1)=25.1;%真实温度初值
P(1)=0.1;%初始协方差
Z(1)=24.9;%初始温度计测量值
Xkf(1)=Z(1);%kalman滤波器初始估计值

%噪声
Q=0.01;
R=0.25;
W=sqrt(Q)*randn(1,N);
V=sqrt(R)*randn(1,N);

%系统矩阵
F=1;
G=1;
H=1;
I=eye(1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%模拟房间温度与测量过程
for k=2:N
    X(k)=F*X(k-1)+G*W(k-1);%k时刻真实温度
    Z(k)=H*X(k)+V(k);%k时刻测量温度
    if Z(k)<=24.5
        Z(k)=24.5;
    end
    X_pre=F*Xkf(k-1);%状态预测           
    P_pre=F*P(k-1)*F'+Q; %协方差预测       
    Kg=P_pre*H'*inv(H*P_pre*H'+R); %kalman增益
    e=Z(k)-H*X_pre;%新息            
    Xkf(k)=X_pre+Kg*e;%状态更新         
    P(k)=(I-Kg*H)*P_pre;%协方差更新
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Err_Messure=zeros(1,N);
Err_Kalman=zeros(1,N);
for k=1:N
    Err_Messure(k)=abs(Z(k)-X(k));%测量与真实值之间的偏差
    Err_Kalman(k)=abs(Xkf(k)-X(k));%估计与真实值之间的偏差
end
t=1:N;
figure('Name','Kalman Filter Simulation','NumberTitle','off');
plot(t,Xexpect,'-b',t,X,'-r',t,Z,'-k',t,Xkf,'-g');
% axis([1,N,24,40]);
legend('expected','real','measure','kalman extimate');         
xlabel('sample time');
ylabel('temperature');
title('Kalman Filter Simulation');
figure('Name','Error Analysis','NumberTitle','off');
plot(t,Err_Messure,'-b',t,Err_Kalman,'-k');
legend('messure error','kalman error');         
xlabel('sample time');
