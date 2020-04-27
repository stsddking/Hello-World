clc;clear

%fn=['D:\S1\1\NS.txt'];  %Point to the path of data
fid=fopen(fn,'r');
[RecorderTime,Nied,EW]=read_remos(fid);

x=EW(13000:14023); 
%xb=[1:12:12000];
%x=x(xb);
N=length(x);
K=20;
ox=x';
ff=abs(fft(ox,N));

% dff=abs(fft(dx,N));
% 
% figure(1);
% subplot(2,1,1)
% hold on;
% plot(ox,'b')
% hold on;
% subplot(2,1,1)
% plot(dx,'r')  
% legend('Original','Noised')
% 
% title('Original and Noised')
% subplot(2,1,2)
% hold on;
% plot(ff,'b')    
% hold on;
% subplot(2,1,2)
% plot(dff,'r')    
% title('Fourier Transform')

snr=0;
mse=0;

snrbut=[];
msebut=[];
snrdct=[];
msedct=[];
snrdft=[];
msedft=[];

elapsed_time = zeros(10,1) ;  % initilaize the elapsed times
fx=[];

 WW1=dctmtx(1024);
 WW2=dctmtx(1024);
 WW3=dctmtx(1024);
 WW4=dctmtx(1024);
 WW=[WW1,WW2,WW3,WW4];
DDCT=sqrt(1/1024)*WW;
for i = 1:4096
    DDCT(:,i) = DDCT(:,i) / norm(DDCT(:,i));
end

WW1=dftmtx(1024);
WW2=dftmtx(1024);
WW3=dftmtx(1024);
WW4=dftmtx(1024);
 WW=[WW1,WW2,WW3,WW4];
DDFT=sqrt(1/1024)*WW;
for i = 1:4096
    DDFT(:,i) = DDFT(:,i) / norm(DDFT(:,i));
end


for avnum=1:10;
dy=randn(N,1);
dy=dy-mean(dy);
dy=dy/std(dy);
dy=dy*sqrt(2^2); %sigma=2 or 10

% dy=dy*sqrt(avnum^2);  %figure 424

dx=ox+dy;

fx=[;matfilter(dx(1:256));matfilter(dx(257:512));matfilter(dx(513:768));matfilter(dx(769:1024))];
snr=10*log10(norm(ox)^2/norm(fx-ox)^2);
mse=(norm(fx-ox))^2/N;
snrbut=[snrbut,snr];
msebut=[msebut,mse];

[AA,RR,indd]=OMP(DDCT,dx,K);
%hat_x=DD*AA;
hat_x=real(DDCT*AA);
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
snrdct=[snrdct,snr];
msedct=[msedct,mse];

[AA,RR,indd]=OMP(DDFT,dx,K);
%hat_x=DD*AA;
hat_x=real(DDFT*AA);
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
snrdft=[snrdft,snr];
msedft=[msedft,mse];

end

figure(1);
hold on;
%plot(hat_x,'k.-','LineWidth',2)                                 
plot(snrbut,'-ro','LineWidth',2);
plot(snrdct,'-go','LineWidth',2);
plot(snrdft,'-bo','LineWidth',2);

legend('FIR-FILTER','DCT','DFT');

figure(2);
hold on;
%plot(hat_x,'k.-','LineWidth',2)                                
plot(msebut,'-ro','LineWidth',2);
plot(msedct,'-go','LineWidth',2);
plot(msedft,'-bo','LineWidth',2);

legend('FIR-FILTER','DCT','DFT');