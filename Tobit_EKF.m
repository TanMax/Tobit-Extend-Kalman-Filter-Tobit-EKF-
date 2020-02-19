%%%TCD 2020/2/7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  函数功能：标量非线性系统扩展Kalman滤波问题
%  状态函数：X(k+1)=X(k) +w(k)
%  观测方程：Z（k）=X(k)^2/10 +v(k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all;
T=75;%采样点
CON=25;%室内温度理论值
Xexpect=CON*ones(1,T);%期望温度
%噪声
Q=0.01;%状态噪声
R=0.01;%观测噪声
w=sqrt(Q)*randn(1,T);
v=sqrt(R)*randn(1,T);
CUT=65;%观测截断值

%初始化
x=zeros(1,T);
x(1)=26;%初始化真实状态值
y=zeros(1,T);
y(1)=x(1)^2/10+v(1);%初始化观测值
P0=0.01;%初始化协方差
Xekf=zeros(1,T);
Xekf(1)=x(1);%初始化kalmam输出初始值

%%%%模拟真实值和观测值
for k=2:T
    x(k)=x(k-1)+w(k-1);
    y(k)=x(k)^2/10+v(k);
    if y(k)<=CUT
        y(k)=CUT;
    end
end

%标准ekf滤波
for k=2:T
    Xn=1*Xekf(k-1);%预测值
    F=1;
    P=F*P0*F'+Q;%%协方差预测
    
    Zn=Xn^2/10;%测量值预测
    H=Xn/5;%观测系数
    
    K=P*H'*inv(H*P*H'+R);%kalman增益
    Xekf(k)=Xn+K*(y(k)-Zn);%状态更新
    P0=(eye(1)-K*H)*P;%%协方差矩阵更新
end

%Tobit ekf滤波
P0_T=P0;
Xekf_T=zeros(1,T);
Xekf_T(1)=x(1)+1;%初始化Tobit ekf
for k=2:T
    Xn_T=1*Xekf_T(k-1);%预测值
    F=1;
    P_T=F*P0_T*F'+Q;%%协方差预测
    
    Zn_T=Xn_T^2/10;%测量值预测
    H_T=Xn_T/5;%观测系数
    q=(Zn_T-CUT)/sqrt(R);%中间变量
    int_pr=cdf('normal',(Zn_T-CUT)/sqrt(R),0,1);%大于截断值的概率
    int_Pr=1-int_pr;%小于截断值的概率
    Pr= normpdf((Zn_T-CUT)/sqrt(R),0,1);%点的概率值
    lamada=Pr/(int_pr);
    Zn1_T=int_pr*(Zn_T+sqrt(R)*lamada)+int_Pr*CUT;%测量值预测更新
    %求kalman增益
    E_Pr=int_pr;%E(p)
    E_Pv=R*(1-lamada*(lamada+q));%E(pvvp)
    Rxy=P_T*H_T'*E_Pr;Ryy=E_Pr*H_T*P_T*H_T'*E_Pr+E_Pv;
    K_T=Rxy*inv(Ryy);
    %%%
    Xekf_T(k)=Xn_T+K_T*(y(k)-Zn1_T);%状态更新
    P0_T=(eye(1)-E_Pr*K_T*H_T)*P_T;%%协方差矩阵更新
end

Xstd=zeros(1,T);
Xstd_T=zeros(1,T);
for k=1:T
    Xstd(k)=abs( Xekf(k)-x(k) );%EKF真实值与估计值偏差
    Xstd_T(k)=abs( Xekf_T(k)-x(k) );%Tobit EKF真实值与估计值偏差
end
t=1:T;
figure('Name','Extend Tobit EKF Filter Simulation','NumberTitle','off');
plot(t,Xexpect,'-b',t,x,'-r',t,Xekf,'-g',t,Xekf_T,'k');
% axis([1,N,24,40]);
legend('expected','real','EKF extimate','Tobit EKF');         
xlabel('sample time');
ylabel('temperature');
title('Extend Tobit Kalman Filter Simulation');
figure('Name','measure','NumberTitle','off');
scatter(t,y);
legend('measure');
