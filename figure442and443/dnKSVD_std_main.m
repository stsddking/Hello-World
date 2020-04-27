clc;clear

%fn=['D:\S1\1\NS.txt'];  %Point to the path of data

fid=fopen(fn,'r');
[RecorderTime,Nied,EW]=read_remos(fid);

x=EW(13000:14023); 
%xb=[1:12:12000];
%x=x(xb);
sigma=2;

N=length(x);
M=fix(N*0.6);
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

snrfir=[];
msefir=[];
snrfj=[];
msefj=[];
snrcs=[];
msecs=[];
snrksvd=[];
mseksvd=[];
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


for avnum=1:10;
    sigma=avnum;
dy=randn(N,1);
dy=dy-mean(dy);
dy=dy/std(dy);
%dy=dy*sqrt(avnum^2);
dy=dy*sigma;
dx=ox+dy;

fx=[;matfilter(dx(1:256));matfilter(dx(257:512));matfilter(dx(513:768));matfilter(dx(769:1024))];
snr=10*log10(norm(ox)^2/norm(fx-ox)^2);
mse=(norm(fx-ox))^2/N;
snrfir=[snrfir,snr];
msefir=[msefir,mse];

[AA,RR,indd]=OMP(DDCT,dx,K);
%hat_x=DD*AA;
hat_x=real(DDCT*AA);
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
snrfj=[snrfj,snr];
msefj=[msefj,mse];

Psi=dctmtx(N);
Phi = randn(M,N);
Phi = orth(Phi')';
A = Phi * Psi;
y=Phi*dx;
tic  
lamda = sigma*sqrt(2*log(N));
theta =  BPDN_quadprog(y,A,lamda);
hat_x = real(Psi * theta);% x=Psi * theta  
toc  
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
snrcs=[snrcs,snr];
msecs=[msecs,mse];

fildNameForGlobalDictionary = 'KSVDTRAIN2';
givenDictionaryFlag = 0;
if (~givenDictionaryFlag)
    eval(['load ',fildNameForGlobalDictionary]);
    DD = Dictionary;
end

[AA,RR,indd]=OMP(DD,dx,K);
hat_x=DD*AA;
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
snrksvd=[snrksvd,snr];
mseksvd=[mseksvd,mse];

end

figure(1);
hold on;
%plot(hat_x,'k.-','LineWidth',2)                                 
plot(snrfir,'-ro','LineWidth',2);
plot(snrfj,'-go','LineWidth',2);
plot(snrcs,'-bo','LineWidth',2);
plot(snrksvd,'-co','LineWidth',2);
%plot(snrdwdft,'-mo','LineWidth',2);
legend('FIR-FILTER','SR','CS','DL');

figure(2);
hold on;
%plot(hat_x,'k.-','LineWidth',2)                                 
plot(msefir,'-ro','LineWidth',2);
plot(msefj,'-go','LineWidth',2);
plot(msecs,'-bo','LineWidth',2);
plot(mseksvd,'-co','LineWidth',2);
%plot(msedwdft,'-mo','LineWidth',2);
legend('FIR-FILTER','SR','CS','DL');