
clc;
clear;

SRCname={'������','��������','ǿ��','����','����','�г�����','��ɸ','�ȶ�','���Ҳ�','����','����','����','����','����'};
load('SRCtrain.mat');
IDtrain=trainmark;

K = 14;   %14��
a=1005; 
aat=[4,4,40,5,10,48,6,6,12,40,40,40,40,40];
aatc=cumsum(aat*3);


A=traindata;
m=600; %ÿ�г���600

%A = double(A);
% change the dimension of d
%d=600;
%R=dctmtx(600);
%R = randn(d,m);
%ORTH_R = orth(R')';
%ORTH_R = ORTH_R';
%RA = ORTH_R*A;

RA=double(A*1);
Coef_Mat = zeros(a,K);   %a��ÿ��Ŀ¼�µ����ļ�����K=61��

RA=sqrt(1/600)*RA;
for i = 1:1005
    RA(:,i) = RA(:,i) / norm(RA(:,i));
end

load('SRCtest.mat');
IDtest=testmark;
T=testdata;
vt=T(:,104);

figure(1);
hold on;
plot(vt,'b');
%y = ORTH_R*vt;
y=double(vt*1);

%K=20;
[AA,RR,indd]=OMP(RA,y,180);
coef=AA;
hat_y=real(RA*AA);
hold on;
plot(hat_y,'r');

for j=1:K
    if j==1
       Coef_Mat(1:aatc(j),j) = coef(1:aatc(j));
    else
        Coef_Mat(aatc(j-1)+1:aatc(j),j) = coef(aatc(j-1)+1:aatc(j));
    end    
end
for j=1:K
    res(j) = norm(y-RA*Coef_Mat(:,j),2);
end
[C,I] = min(res);

title(['The signal type is   -----  ','number��',num2str(I),'    ',SRCname(I)]);

            

    
    