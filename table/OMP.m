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
residual=X; %��ʼ���в�
indx=zeros(L,1);
for i=1:L
    proj=D'*residual;%Dת����residual��ˣ��õ���residual��Dÿһ�е��ڻ�ֵ
    [~,pos]=max(abs(proj));%�ҵ��ڻ����ֵ��λ��
    pos=pos(1);%�����ֵ��ֹһ����ȡ��һ��
    indx(i)=pos;%�����λ�ô����������ĵ�j��ֵ
    a=pinv(D(:,indx(1:i)))*X;%indx(1:j)��ʾ��һ��ǰj��Ԫ��
    residual=X-D(:,indx(1:i))*a;
    res=norm(residual);
    if i==1
    threshold=res;
    end
    if res< threshold*0.001
        break;
    end
end
A=zeros(size(D,2),1);
A(indx(indx~=0))=a;
ind=indx;
end