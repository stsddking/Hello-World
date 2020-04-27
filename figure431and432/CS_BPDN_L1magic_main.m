clc;clear

%fn=['D:\S1\1\NS.txt'];  %Point to the path of data
fid=fopen(fn,'r');
[RecorderTime,Nied,EW]=read_remos(fid);

%%  1. ʱ������ź�����
x=EW(12000:12255); %֡��2�룬�ָ�Ч�����
%xb=[1:12:12000];
%x=x(xb);
N=length(x);
M=fix(N*0.6);
ox=x';
%x=matfilter(ox);
ff=abs(fft(ox,N));
%dx=ox+(randn(N,1)-0.5)*2;
dx=ox+randn(N,1)*0.6;
dff=abs(fft(dx,N));

figure(1);
subplot(2,1,1)
hold on;
plot(ox,'b')
hold on;
subplot(2,1,1)
plot(dx,'r')  
%hold on;
%plot(x,'r')  
legend('Original','Noised')
%title('Original and Filtered')
title('Original and Noised')
subplot(2,1,2)
hold on;
plot(ff,'b')    
hold on;
subplot(2,1,2)
plot(dff,'r')    
title('Fourier Transform')


%figure(2);

%hold on;
%plot(x,'r')  
%legend('Original','Filtered')
%title('Original and Filtered')
%title('Noised')

%title('Fourier Transform')


%[mk,indexff]=max(ff);
%K=sum(ff>mk*0.1);      %  ϡ���(��FFT���Կ�����)
%K=fix(N*0.1);
K=30;
%M=round(K*log(N/K));     %  ������(M>=K*log(N/K),����40,���г���ĸ���)
%Psi = eye(N);%x������ϡ��ģ�����ϡ�����Ϊ��λ��x=Psi*theta  
Psi=dctmtx(N);
Phi = randn(M,N);%��������Ϊ��˹����
Phi = orth(Phi')';
A = Phi * Psi;%���о���  

%y = Phi * ox + e;%�õ��۲�����y  
y=Phi*dx;
x0 = A'*y;
%% �ָ��ع��ź�x  
sigma = 0.6;
%e = sigma*randn(M,1);
epsilon =  sigma*sqrt(M)*sqrt(1 + 2*sqrt(2)/sqrt(M));
tic  
%lamda = sigma*sqrt(2*log(N));
%theta =  BPDN_quadprog(y,A,lamda);
theta= l1qc_logbarrier(x0 , A, [], y, epsilon, 1e-3);
hat_x = Psi * theta;% x=Psi * theta  
toc  



figure(2);
hold on;
plot(hat_x,'k.-','LineWidth',2)                                 %  �ؽ��ź�
hold on;
plot(ox,'r')                                       %  ԭʼ�ź�
%hold on;
%plot(dx,'g')                                       %  �����ź�
%legend('Recovery','Original','Noised')
%norm(hat_x.'-x)/norm(x)                           %  �ع����
%fprintf('\n�ָ��в');  
%tt=norm(hat_x-x.');%�ָ��в�
%title(['�ָ��в',num2str(tt)]);
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
%erp=(norm(hat_x)^2/norm(x)^2)*100;
%title(['SNR��',num2str(snr),'  ','MSE(%)��',num2str(mse),'  ','ERP(%)��',num2str(erp)]);
title(['BP Denoising   -----  ','SNR��',num2str(snr),'    ','MSE��',num2str(mse)]);

fx=[];
%fx=[;matfilter(dx(1:256));matfilter(dx(257:512));matfilter(dx(513:768));matfilter(dx(769:1024))];
fx=matfilter(dx);
figure(3);
hold on;
plot(fx,'k.-','LineWidth',2)                                 %  �ؽ��ź�
hold on;
plot(ox,'r')                                       %  ԭʼ�ź�
hold on;
plot(dx,'g')                                       %  �����ź�
legend('filtered','Original','Noised')
%norm(hat_x.'-x)/norm(x)                           %  �ع����
%fprintf('\n�ָ��в');  
%tt=norm(hat_x-x.');%�ָ��в�
%title(['�ָ��в',num2str(tt)]);
snr=10*log10(norm(ox)^2/norm(fx-ox)^2);
mse=(norm(fx-ox))^2/N;
%erp=(norm(hat_x)^2/norm(x)^2)*100;
%title(['SNR��',num2str(snr),'  ','MSE(%)��',num2str(mse),'  ','ERP(%)��',num2str(erp)]);
title(['Filered   -----  ','SNR��',num2str(snr),'    ','MSE��',num2str(mse)]);