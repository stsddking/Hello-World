clc;clear

%fn=['D:\S1\1\NS.txt'];  %Point to the path of data

fid=fopen(fn,'r');
[RecorderTime,Nied,EW]=read_remos(fid);
sigma=10;
%-------------------------------------------------------------------------------------------------------------------

%%  1. ʱ������ź�����------------------------------------------------------------------------------------
x=EW(13000:14023); %֡��2�룬�ָ�Ч�����
%xb=[1:12:12000];
%x=x(xb);
N=length(x);

ox=x';
%x=matfilter(ox);
ff=abs(fft(ox,N));
%dx=ox+(rand(N,1)-0.5)*4;
%dx=awgn(ox,-10,'measured');
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


%------------------------------------------------------------------------------------------------------------------

K=20;
%M=round(K*log(N/K));     %  ������(M>=K*log(N/K),����40,���г���ĸ���)
%DD=OCT_DCTdic(16,32);
fildNameForGlobalDictionary = 'KSVDTRAIN2';
givenDictionaryFlag = 0;
if (~givenDictionaryFlag)
    eval(['load ',fildNameForGlobalDictionary]);
    DD = Dictionary;
end

[AA,RR,indd]=OMP(DD,dx,K);
hat_x=DD*AA;

figure(2);
hold on;
plot(hat_x,'k.-','LineWidth',2)                                 %  �ؽ��ź�
hold on;
plot(ox,'r')                                       %  ԭʼ�ź�
hold on;
plot(dx,'g')                                       %  �����ź�
legend('Recovery','Original','Noised')
%norm(hat_x.'-x)/norm(x)                           %  �ع����
%fprintf('\n�ָ��в');  
%tt=norm(hat_x-x.');%�ָ��в�
%title(['�ָ��в',num2str(tt)]);
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
%erp=(norm(hat_x)^2/norm(x)^2)*100;
%title(['SNR��',num2str(snr),'  ','MSE(%)��',num2str(mse),'  ','ERP(%)��',num2str(erp)]);
title(['OMP Denoising   -----  ','SNR��',num2str(snr),'    ','MSE��',num2str(mse)]);

fx=[];
fx=[;matfilter(dx(1:256));matfilter(dx(257:512));matfilter(dx(513:768));matfilter(dx(769:1024))];
%fx=matfilter(dx);
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
title(['Filtered   -----  ','SNR��',num2str(snr),'    ','MSE��',num2str(mse)]);