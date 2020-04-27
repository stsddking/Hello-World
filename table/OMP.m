%OMP计算稀疏系数
function [A,res,ind]=OMP(D,X,L)
% 输入参数:
%       D - 过完备字典，注意：必须字典的各列必须经过了规范化
%       X - 信号
%       L - 稀疏度，系数中非零元个数的最大值
% 输出参数:
%       A - 当前信号的系数
%       res - 残差
 
%%
residual=X; %初始化残差
indx=zeros(L,1);
for i=1:L
    proj=D'*residual;%D转置与residual相乘，得到与residual与D每一列的内积值
    [~,pos]=max(abs(proj));%找到内积最大值的位置
    pos=pos(1);%若最大值不止一个，取第一个
    indx(i)=pos;%将这个位置存入索引集的第j个值
    a=pinv(D(:,indx(1:i)))*X;%indx(1:j)表示第一列前j个元素
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