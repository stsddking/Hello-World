

clc;
clear;

SRCname={'白噪声','背景噪声','强夯','捣固','方波','列车干扰','清筛','稳定','正弦波','地震波','地震波','地震波','地震波','地震波'};
load('SRCtrain.mat');
IDtrain=trainmark;
SRCresult=[];
K = 14;   %14类
a=1005; 
aat=[4,4,40,5,10,48,6,6,12,40,40,40,40,40];
aatc=cumsum(aat*3);
bbt=[1,1,10,1,3,12,1,1,4,10,10,10,10,10];
bbtc=cumsum(bbt*3);



A=traindata;
m=600; %每列长度600

%A = double(A);
% change the dimension of d
%d=600;
%R = randn(d,m);
%ORTH_R = orth(R');
%ORTH_R = ORTH_R';

%RA = ORTH_R*A;
RA=A;
Coef_Mat = zeros(a,K);   %a是每个目录下的总文件数，K=61？

 RA=sqrt(1/600)*RA;
 for i = 1:1005
     RA(:,i) = RA(:,i) / norm(RA(:,i));
 end

load('SRCtest.mat');
IDtest=testmark;
T=testdata;
for ii=1:252
vt=T(:,ii);

% figure(1);
% hold on;
% plot(vt,'b');
%y = ORTH_R*vt;
y=vt;
    

[AA,RR,indd]=OMP(RA,y,200);
coef=AA;
hat_y=real(RA*AA);
% hold on;
% plot(hat_y,'r');

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
SRCresult(ii)=I;
end
% title(['The signal type is   -----  ','number：',num2str(I),'    ',SRCname(I)]);
rwj=[];
 ccc=IDtest(1,:);
for ii=1:252
    if SRCresult(ii)==ccc(ii);
        rwj(ii)=1;
    else
        rwj(ii)=0;
    end   
end    

rwjz=rwj(1:102);
rwjzt=find(rwjz(:)~=0);
zzr=length(rwjzt)/length(rwjz);
rjz=SRCresult(1:102);
rwjzt2=find(rjz(:)<10);
zzr2=length(rwjzt2)/length(rjz);

rwjd=rwj(103:252);
rwjdt=find(rwjd(:)~=0);
dzr=length(rwjdt)/length(rwjd);
rjd=SRCresult(103:252);
rwjdt2=find(rjd(:)>9);
dzr2=length(rwjdt2)/length(rjd);