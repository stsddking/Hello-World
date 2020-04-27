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
K=20;
sigma = 1;

snr=0;
mse=0;

snrfir=[];
msefir=[];
snrsromp=[];
msesromp=[];
snrsrbpdn=[];
msesrbpdn=[];
snrcsomp=[];
msecsomp=[];
snrcsbpdn=[];
msecsbpdn=[];

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
WW1=dctmtx(1024);
WW2=dctmtx(1024);
WW3=dctmtx(1024);
WW4=dctmtx(1024);
WW=[WW1,WW2,WW3,WW4]; 
DD=sqrt(1/1024)*WW;
for i = 1:4096
    DD(:,i) = DD(:,i) / norm(DD(:,i));
end

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


[AA,RR,indd]=OMP(DD,dx,K);
hat_x=real(DD*AA);
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
snrsromp=[snrsromp,snr];
msesromp=[msesromp,mse];


%DD = orth(DD')';
tic  
lamda = sigma*sqrt(2*log(N));
theta =  BPDN_quadprog(dx,DD,lamda);
hat_x = DD* theta;% x=Psi * theta  
toc  
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
snrsrbpdn=[snrsrbpdn,snr];
msesrbpdn=[msesrbpdn,mse];

Psi=dctmtx(N);
Phi = randn(M,N);%测量矩阵为高斯矩阵
Phi = orth(Phi')';
s=Phi*dx;
%% 恢复重构信号x  
tic  
T=Phi*Psi';                                       %  恢复矩阵(测量矩阵*正交反变换矩阵)
hat_y=zeros(1,N);                                 %  待重构的谱域(变换域)向量                     
Aug_t=[];                                         %  增量矩阵(初始值为空矩阵)
r_n=s;                                            %  残差值
for times=1:K;                                    %  迭代次数(有噪声的情况下,该迭代次数为K)
    for col=1:N;                                  %  恢复矩阵的所有列向量
        product(col)=abs(T(:,col)'*r_n);          %  恢复矩阵的列向量和残差的投影系数(内积值) 
    end
    [val,pos]=max(product);                       %  最大投影系数对应的位置
    Aug_t=[Aug_t,T(:,pos)];                       %  矩阵扩充
    T(:,pos)=zeros(M,1);                          %  选中的列置零（实质上应该去掉，为了简单我把它置零）
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s;           %  最小二乘,使残差最小
    r_n=s-Aug_t*aug_y;                            %  残差
    pos_array(times)=pos;                         %  纪录最大投影系数的位置
end
hat_y(pos_array)=aug_y;                           %  重构的谱域向量
hat_x=real(Psi'*hat_y.');                         %  做逆傅里叶变换重构得到时域信号
toc  
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
snrcsomp=[snrcsomp,snr];
msecsomp=[msecsomp,mse];


Psi=dctmtx(N);
Phi = randn(M,N);%测量矩阵为高斯矩阵
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
snrcsbpdn=[snrcsbpdn,snr];
msecsbpdn=[msecsbpdn,mse];

end

figure(1);
hold on;
%plot(hat_x,'k.-','LineWidth',2)                                 %  重建信号
plot(snrfir,'-ro','LineWidth',2);
plot(snrsromp,'-go','LineWidth',2);
plot(snrsrbpdn,'-bo','LineWidth',2);
plot(snrcsomp,'-co','LineWidth',2);
plot(snrcsbpdn,'-mo','LineWidth',2);
%plot(snrsr,'-yo','LineWidth',2);
legend('FIR-FILTER','SR-OMP','SR-BPDN','CS-OMP','CS-BPDN');


figure(2);
hold on;
plot(msefir,'-ro','LineWidth',2);
plot(msesromp,'-go','LineWidth',2);
plot(msesrbpdn,'-bo','LineWidth',2);
plot(msecsomp,'-co','LineWidth',2);
plot(msecsbpdn,'-mo','LineWidth',2);
%plot(msesr,'-yo','LineWidth',2);%plot(hat_x,'k.-','LineWidth',2)                                 %  重建信号
legend('FIR-FILTER','SR-OMP','SR-BPDN','CS-OMP','CS-BPDN');

