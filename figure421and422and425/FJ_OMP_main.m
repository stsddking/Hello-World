clc;clear

%fn=['D:\S1\1\NS.txt'];  %Point to the path of data
fid=fopen(fn,'r');
[RecorderTime,Nied,EW]=read_remos(fid);


x=EW(13000:14023); 

N=length(x);

ox=x';

ff=abs(fft(ox,N));

sigma = 10;
dy=randn(N,1);
dy=dy-mean(dy);
dy=dy/std(dy);
dy=dy*sigma;
dx=ox+dy;

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
%K=sum(ff>mk*0.1);      %  稀疏度(做FFT可以看出来)
%K=fix(N*0.1);
K=15;
%M=round(K*log(N/K));     %  测量数(M>=K*log(N/K),至少40,但有出错的概率)
%DD=OCT_DCTdic(32,64);
%DD=sqrt(1/1024)*OCT_DCTdic(32,64);

%WW=sqrt(1/1024)*randn(1024,4096);
WW1=dctmtx(1024);
WW2=dctmtx(1024);
WW3=dctmtx(1024);
WW4=dctmtx(1024);

% WW1=dftmtx(1024);
% WW2=dftmtx(1024);
% WW3=dftmtx(1024);
% WW4=dftmtx(1024);

WW=[WW1,WW2,WW3,WW4]; 

DD=sqrt(1/1024)*WW;
for i = 1:4096
    DD(:,i) = DD(:,i) / norm(DD(:,i));
end
[AA,RR,indd]=OMP(DD,dx,K);
%hat_x=DD*AA;
hat_x=real(DD*AA);

figure(2);
hold on;
plot(hat_x,'k.-','LineWidth',2)                                 %  重建信号
hold on;
plot(ox,'r')                                       %  原始信号
hold on;
plot(dx,'g')                                       %  噪声信号
legend('Recovery','Original','Noised')
%norm(hat_x.'-x)/norm(x)                           %  重构误差
%fprintf('\n恢复残差：');  
%tt=norm(hat_x-x.');%恢复残差
%title(['恢复残差：',num2str(tt)]);
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
%erp=(norm(hat_x)^2/norm(x)^2)*100;
%title(['SNR：',num2str(snr),'  ','MSE(%)：',num2str(mse),'  ','ERP(%)：',num2str(erp)]);
title(['OMP Denoising   -----  ','SNR：',num2str(snr),'    ','MSE：',num2str(mse)]);

fx=[];
fx=[;matfilter(dx(1:256));matfilter(dx(257:512));matfilter(dx(513:768));matfilter(dx(769:1024))];
%=matfilter(dx);
figure(3);
hold on;
plot(fx,'k.-','LineWidth',2)                                 %  重建信号
hold on;
plot(ox,'r')                                       %  原始信号
hold on;
plot(dx,'g')                                       %  噪声信号
legend('filtered','Original','Noised')
%norm(hat_x.'-x)/norm(x)                           %  重构误差
%fprintf('\n恢复残差：');  
%tt=norm(hat_x-x.');%恢复残差
%title(['恢复残差：',num2str(tt)]);
snr=10*log10(norm(ox)^2/norm(fx-ox)^2);
mse=(norm(fx-ox))^2/N;
%erp=(norm(hat_x)^2/norm(x)^2)*100;
%title(['SNR：',num2str(snr),'  ','MSE(%)：',num2str(mse),'  ','ERP(%)：',num2str(erp)]);
title(['Filtered   -----  ','SNR：',num2str(snr),'    ','MSE：',num2str(mse)]);