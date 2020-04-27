clc;clear

%fn=['D:\S1\1\NS.txt'];  %Point to the path of data

fid=fopen(fn,'r');
[RecorderTime,Nied,EW]=read_remos(fid);

%%  1. ʱ������ź�����
x=EW(13000:14023); %֡��2�룬�ָ�Ч�����
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
Phi = randn(M,N);%��������Ϊ��˹����
Phi = orth(Phi')';
s=Phi*dx;
%% �ָ��ع��ź�x  
tic  
T=Phi*Psi';                                       %  �ָ�����(��������*�������任����)
hat_y=zeros(1,N);                                 %  ���ع�������(�任��)����                     
Aug_t=[];                                         %  ��������(��ʼֵΪ�վ���)
r_n=s;                                            %  �в�ֵ
for times=1:K;                                    %  ��������(�������������,�õ�������ΪK)
    for col=1:N;                                  %  �ָ����������������
        product(col)=abs(T(:,col)'*r_n);          %  �ָ�������������Ͳв��ͶӰϵ��(�ڻ�ֵ) 
    end
    [val,pos]=max(product);                       %  ���ͶӰϵ����Ӧ��λ��
    Aug_t=[Aug_t,T(:,pos)];                       %  ��������
    T(:,pos)=zeros(M,1);                          %  ѡ�е������㣨ʵ����Ӧ��ȥ����Ϊ�˼��Ұ������㣩
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*s;           %  ��С����,ʹ�в���С
    r_n=s-Aug_t*aug_y;                            %  �в�
    pos_array(times)=pos;                         %  ��¼���ͶӰϵ����λ��
end
hat_y(pos_array)=aug_y;                           %  �ع�����������
hat_x=real(Psi'*hat_y.');                         %  ���渵��Ҷ�任�ع��õ�ʱ���ź�
toc  
snr=10*log10(norm(ox)^2/norm(hat_x-ox)^2);
mse=(norm(hat_x-ox))^2/N;
snrcsomp=[snrcsomp,snr];
msecsomp=[msecsomp,mse];


Psi=dctmtx(N);
Phi = randn(M,N);%��������Ϊ��˹����
Phi = orth(Phi')';
A = Phi * Psi;%���о���  
y=Phi*dx;
%% �ָ��ع��ź�x  
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
%plot(hat_x,'k.-','LineWidth',2)                                 %  �ؽ��ź�
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
%plot(msesr,'-yo','LineWidth',2);%plot(hat_x,'k.-','LineWidth',2)                                 %  �ؽ��ź�
legend('FIR-FILTER','SR-OMP','SR-BPDN','CS-OMP','CS-BPDN');

