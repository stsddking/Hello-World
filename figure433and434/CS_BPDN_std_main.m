clc;clear
%fn=['D:\S1\1\NS.txt'];  %Point to the path of data

fid=fopen(fn,'r');
[RecorderTime,Nied,EW]=read_remos(fid);

%%  1. 时域测试信号生成
x=EW(13000:14023); %帧长2秒，恢复效果最好
%xb=[1:12:12000];
%x=x(xb);
N=length(x);
M=fix(N*0.6);
ox=x';
%x=matfilter(ox);
ff=abs(fft(ox,N));
%dx=ox+(randn(N,1)-0.5)*2;
%dx=ox+rand(N,1)*2;
fx=[];
K=30;
sigma = 1;

snr=0;
mse=0;

snrfir=[];
msefir=[];
snrgs=[];
msegs=[];
snrbnl=[];
msebnl=[];
snrhdm=[];
msehdm=[];
snrtpl=[];
msetpl=[];
snrsr=[];
msesr=[];
% figure(1);
% subplot(2,1,1)
% hold on;
% plot(ox,'b')
% hold on;
% subplot(2,1,1)
% plot(dx,'r')  
% legend('Original','Noised')
% title('Original and Noised')
% subplot(2,1,2)
% hold on;
% plot(ff,'b')    
% hold on;
% subplot(2,1,2)
% plot(dff,'r')    
% title('Fourier Transform')
for avnum=1:10;
sigma=avnum;
dy=randn(N,1);
dy=dy-mean(dy);
dy=dy/std(dy);
dy=dy*sigma;
dx=ox+dy;
fx=[;matfilter(dx(1:256));matfilter(dx(257:512));matfilter(dx(513:768));matfilter(dx(769:1024))];
snr=10*log10(norm(ox)^2/norm(fx-ox)^2);
mse=(norm(fx-ox))^2/N;
snrfir=[snrfir,snr];
msefir=[msefir,mse];

%dff=abs(fft(dx,N));
%Psi=real(dftmtx(N));
Psi=dctmtx(N);
Phi = randn(M,N);%测量矩阵为高斯矩阵
%Phi=BernoulliMtx(M,N);
%Phi=PartFourierMtx(M,N);
%Phi=PartHadamardMtx(M,N);
%Phi=ToeplitzMtx(M,N);
%Phi=SparseRandomMtx(M,N,150);
Phi = orth(Phi')';
A = Phi * Psi;%传感矩阵  
%sigma = 2;
%e = sigma*randn(M,1);
%y = Phi * ox + e;%得到观测向量y  
y=Phi*dx;
%y=A*dx;
%% 恢复重构信号x  
tic  
lamda = sigma*sqrt(2*log(N));
theta =  BPDN_quadprog(y,A,lamda);
hat_x = real(Psi * theta);% x=Psi * theta  
toc  
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
snrgs=[snrgs,snr];
msegs=[msegs,mse];

Psi=dctmtx(N);
Phi=BernoulliMtx(M,N);
%Phi=PartFourierMtx(M,N);
%Phi=PartHadamardMtx(M,N);
%Phi=ToeplitzMtx(M,N);
%Phi=SparseRandomMtx(M,N,150);
Phi = orth(Phi')';
A = Phi * Psi;%传感矩阵  
y=Phi*dx;
%% 恢复重构信号x  
tic  
lamda = sigma*sqrt(2*log(N));
theta =  BPDN_quadprog(y,A,lamda);
hat_x = real(Psi * theta);% x=Psi * theta  
toc  
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
snrbnl=[snrbnl,snr];
msebnl=[msebnl,mse];

Psi=dctmtx(N);
Phi=PartHadamardMtx(M,N);
%Phi=ToeplitzMtx(M,N);
%Phi=SparseRandomMtx(M,N,150);
Phi = orth(Phi')';
A = Phi * Psi;%传感矩阵  
y=Phi*dx;
%% 恢复重构信号x  
tic  
lamda = sigma*sqrt(2*log(N));
theta =  BPDN_quadprog(y,A,lamda);
hat_x = real(Psi * theta);% x=Psi * theta  
toc  
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
snrhdm=[snrhdm,snr];
msehdm=[msehdm,mse];

Psi=dctmtx(N);
Phi=ToeplitzMtx(M,N);
%Phi=SparseRandomMtx(M,N,150);
Phi = orth(Phi')';
A = Phi * Psi;%传感矩阵  
y=Phi*dx;
%% 恢复重构信号x  
tic  
lamda = sigma*sqrt(2*log(N));
theta =  BPDN_quadprog(y,A,lamda);
hat_x = real(Psi * theta);% x=Psi * theta  
toc  
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
snrtpl=[snrtpl,snr];
msetpl=[msetpl,mse];

Psi=dctmtx(N);
Phi=SparseRandomMtx(M,N,150);
Phi = orth(Phi')';
A = Phi * Psi;%传感矩阵  
y=Phi*dx;
%% 恢复重构信号x  
tic  
lamda = sigma*sqrt(2*log(N));
theta =  BPDN_quadprog(y,A,lamda);
hat_x = real(Psi * theta);% x=Psi * theta  
toc  
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
snrsr=[snrsr,snr];
msesr=[msesr,mse];
end

figure(1);
hold on;
%plot(hat_x,'k.-','LineWidth',2)                                 %  重建信号
plot(snrfir,'-ro','LineWidth',2);
plot(snrgs,'-go','LineWidth',2);
plot(snrbnl,'-bo','LineWidth',2);
plot(snrhdm,'-co','LineWidth',2);
plot(snrtpl,'-mo','LineWidth',2);
plot(snrsr,'-yo','LineWidth',2);
legend('FIR-FILTER','GuassMtx','BernoulliMtx','PartHadamardMtx','ToeplitzMtx','SparseRandomMtx');


figure(2);
hold on;
plot(msefir,'-ro','LineWidth',2);
plot(msegs,'-go','LineWidth',2);
plot(msebnl,'-bo','LineWidth',2);
plot(msehdm,'-co','LineWidth',2);
plot(msetpl,'-mo','LineWidth',2);
plot(msesr,'-yo','LineWidth',2);%plot(hat_x,'k.-','LineWidth',2)                                 %  重建信号
legend('FIR-FILTER','GuassMtx','BernoulliMtx','PartHadamardMtx','ToeplitzMtx','SparseRandomMtx');
