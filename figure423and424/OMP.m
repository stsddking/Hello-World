%OMP����ϡ��ϵ��
function [A,res,ind]=OMP(D,X,L)
% �������:
%       D - ���걸�ֵ䣬ע�⣺�����ֵ�ĸ��б��뾭���˹淶��
%       X - �ź�
%       L - ϡ��ȣ�ϵ���з���Ԫ���������ֵ
% �������:
%       A - ��ǰ�źŵ�ϵ��
%       res - �в�
 
%%
M=size(D,1);
Aug_t=[];   
residual=X; %��ʼ���в�
indx=zeros(L,1);
for i=1:L
    %proj=D'*residual;%Dת����residual��ˣ��õ���residual��Dÿһ�е��ڻ�ֵ
    jn=size(D,2);
    for j=1:jn
       Dj=[D(:,j)';residual'];
        proj(j)=pdist(Dj,'euclidean');
        %proj(j)=pdist(Dj,'cityblock');
        %proj(j)=pdist(Dj,'chebychev');
        %proj(j)=pdist(Dj,'minkowski',2);
        %proj(j)=pdist(Dj,'seuclidean',[0.5,1]);
        %proj(j)=pdist(Dj,'mahalanobis');
        %proj(j)=1-pdist(Dj,'cosine');
        %proj(j)=pdist(Dj,'correlation');
    end
      % [~,pos]=max(abs(proj));%�ҵ��ڻ����ֵ��λ��
    %[~,pos]=max(proj);
    [~,pos]=min(abs(proj));
    pos=pos(1);%�����ֵ��ֹһ����ȡ��һ��
    indx(i)=pos;%�����λ�ô����������ĵ�j��ֵ
     Aug_t=[Aug_t,D(:,pos)];  
     D(:,pos)=zeros(M,1);
    %a=pinv(D(:,indx(1:i)))*X;%indx(1:j)��ʾ��һ��ǰj��Ԫ��
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*X; 
    residual=X-Aug_t*aug_y;
    %res=norm(residual);
    res=residual;
    if i==1
    threshold=res;
    end
   if res< threshold*0.001
    %if res< 0.001
        break;
    end
end
A=zeros(size(D,2),1);
%A(indx(indx~=0))=a;
A(indx(indx~=0))=aug_y;
ind=indx;
end