clear all;		
close all;
clc;
clf;%用来清除图形的命令

Nrun =10;%独立运行次数

L = 14;    %滤波器阶数为20阶，对应于volterra滤波器的阶数14阶
M = 3000;%迭代次数
lamda=0;  %0.01;
% u = 0.005;  %LMP的步长0.005

NKE1 = zeros(1,M);
NKE2 = zeros(1,M);
% EP = zeros(1,M);
MSD_all = zeros(1,M);
K=4;

Hk=[1.00 -0.80 0 1.90 0.95 0 1.10 0 0 0 0 0 -0.63 0]';%待辨识的二阶非线性系统的相应权值
for j=1:Nrun
%实验条件和Nonlinear system identification in impulsive environments(2005)一样
       
   x = randn(1,M); %高斯白噪声,Ex1产生高斯分布的信号
%    x = unifrnd(-1,1,1,M); %高斯白噪声,Ex2产生均匀分布的信号
   xl1=zeros(1,M);
   xl2=xl1;
   xl3=xl1;
   xl1(2:M)=x(1:M-1);
   xl2(3:M)=x(1:M-2);
   xl3(4:M)=x(1:M-3); 

   d=zeros(1,M);

     d1 = 1*x - 0.8*xl1 + 0*xl2 + 1.9*xl3 + ... %一次项
       0.95*x.^2 + 0*x.*xl1 + 1.1*x.*xl2 + 0*x.*xl3 + ... 
       0*xl1.^2 + 0*xl1.*xl2 + 0*xl1.*xl3 + ...
       0*xl2.^2 - 0.63*xl2.*xl3 + 0*xl3.^2;%二次项
%      d2=awgn(d1,20);%加入信噪比为30dB的高斯白噪声,实验结果表明是否加入高斯白噪声对结果没有影响
     d2 = d1;
     a_stable = rasd(M,1.6,0,1/15,0);
     d = d2 + a_stable;%a = 1.25时能够达到接近-30dB，与2005trans的论文具有相似的结果


   w=zeros(L,M);           %RLogLMS算法权向量
        
   uxl=[x;xl1;xl2;xl3;x.^2;x.*xl1;x.*xl2;x.*xl3;xl1.^2;xl1.*xl2;xl1.*xl3;xl2.^2; ... 
       xl2.*xl3;xl3.^2]; %输入向量
   
   forget_factor = 0.995;
   del = 0.01;%正则化参数
   p = (1/del)*eye(L,L);%卡尔曼增益矩阵
   
   for i=K:M
       %------------------RLogLMP---------------------
      e(i) = d(i)-w(:,i-1)'*uxl(:,i); % error sample
      pp = 0.8;%%
%       vv = (log(abs(e(i))))^(pp-1)/abs(e(i));%原始论文提出的代价函数，由于没有绝对值，只能够适用于整数次方的形式
      vv = abs( (log(abs(e(i))) )^(pp-1))/(abs(e(i)));%论文改4中的形式
%       vv = ( sign(log(abs(e(i)))^pp)*log(abs(e(i)))^(pp-1) ) / (abs(e(i)))^2;
%       vv = abs((log(abs(e(i))))^(pp))/((abs(e(i)))^2) / log(abs(e(i)));
%       vv = abs((log(abs(e(i))))^(pp)) / (((abs(e(i)))^2) * log(abs(e(i))));
      k = vv*p*uxl(:,i) / (forget_factor + vv*uxl(:,i)'*p*uxl(:,i));%%即为z（n）
      w(:,i) = w(:,i-1) + k*e(i);
      p = (1/forget_factor)*(p - k*uxl(:,i)'*p);%%即为P(n)
        
      NKE1(i) = NKE1(i) + norm(w(1:4,i) - [1;-0.8;0;1.9])^2/ norm([1;-0.8;0;1.9])^2;
      NKE2(i) = NKE2(i) + norm((w(5:14,i) - [0.95;0;1.10;0;0;0;0;0;-0.63;0]),'fro')^2/ norm([0.95;0;1.10;0;0;0;0;0;-0.63;0],'fro')^2;
%       EP(i) = EP(i) + e(i);
      MSD_all(i) = MSD_all(i) + 10*log10(norm(w(:,i)-Hk)^2);
   end

%-LMP-
h1(j,:)=w(1,:);
h3(j,:)=w(3,:);
h4(j,:)=w(4,:);
h7(j,:)=w(7,:);
h10(j,:)=w(10,:);
h13(j,:)=w(13,:);
   
end

NKE1 = NKE1/Nrun;
NKE2 = NKE2/Nrun;
MSD_all = MSD_all / Nrun;

%LMP
h1M=mean(h1);
h3M=mean(h3);
h4M=mean(h4);
h7M=mean(h7);
h10M=mean(h10);
h13M=mean(h13);

figure(1);
plot(1:M,MSD_all,'b','LineWidth', 2);
xlabel('\itn','FontSize',12,'FontName','Times New Roman');
ylabel('MSD','FontSize',12,'FontName','Times New Roman');
legend('RLogLMP');

figure(2)
plot(d,'g:','LineWidth', 2);
hold on;
plot(a_stable,'r','LineWidth', 2)
hold off;
xlabel('\itn','FontSize',12,'FontName','Times New Roman');
legend('Desire signal','a-stable(a=1.25)');

figure(3)
semilogy((NKE1),'b','LineWidth', 2);
hold on;
semilogy((NKE2),'r','LineWidth', 2);
hold off;
legend('RLogLMP(h1)','RLogLMP(h2)');
xlabel('\itn','FontSize',12,'FontName','Times New Roman');
ylabel('Normalized Volterra Kernel Error','FontSize',12,'FontName','Times New Roman');

% figure(4)
% xx = [w(1:4,i)',zeros(1,5996)];
% yy = [w(5:14,i)',zeros(1,5990)];
% zz = (abs(EP)).^pp;
% A = [xx;yy;zz];
% xlabel('h1','FontSize',16);
% ylabel('h2','FontSize',16);
% zlabel('J(n)','FontSize',16);
% mesh(A)


% figure(4)
% subplot(2,1,1)
% plot(1:M,h1M(1:M),'b',1:M,1,'g-.');
% xlabel('N','FontSize',16);
% ylabel('h1(n)','FontSize',16);
% % gtext('desired value:1.00'); 
% 
% subplot(2,1,2)
% plot(1:M,h3M(1:M),'b',1:M,0,'g-.');
% xlabel('N','FontSize',16);
% ylabel('h3(n)','FontSize',16);
% % gtext('desired value:0.00'); 
% 
% figure(5)
% subplot(2,1,1)
% plot(1:M,h4M(1:M),'b',1:M,1.9,'g-.');
% xlabel('N','FontSize',16);
% ylabel('h4(n)','FontSize',16);
% % gtext('desired value:1.90'); 
% 
% subplot(2,1,2)
% plot(1:M,h7M(1:M),'b',1:M,1.1,'g-.');
% xlabel('N','FontSize',16);
% ylabel('h7(n)','FontSize',16);
% gtext('desired value:1.10'); 
% save h1_RLogLMP.mat h1M
% save h3_RLogLMP.mat h3M
% save h4_RLogLMP.mat h4M
% save h7_RLogLMP.mat h7M
% save h10_RLogLMP.mat h10M
% save h13_RLogLMP.mat h13M

% save NKE1_RLogLMP.mat NKE1
% save NKE2_RLogLMP.mat NKE2
% save MSD_RLogLMP.mat MSD_all
