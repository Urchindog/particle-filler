% % 随机采样子函数
% function outIndex = randomR(inIndex,q)
% if nargin < 2
%     error('Not enough input arguments.'); 
% end
% outIndex=zeros(size(inIndex));
% [num,col]=size(q);
% u=rand(num,1);
% u=sort(u);
% l=cumsum(q);
% i=1;
% for j=1:num
%     while (i<=num)&(u(i)<=l(j))
%         outIndex(i)=j;
%         i=i+1;
%     end
% end

% function outIndex = randomR(inIndex,q)
% if nargin < 2
%     error('Not enough input arguments.'); 
% end
% outIndex=zeros(size(inIndex));
% [num,col]=size(q);
% u=rand(num,1);
% u=sort(u);
% l=cumsum(q);
%   for i = 1 : N    
%         outIndex(i) = q(find(rand <= cumsum(q),1));   % 粒子权重大的将多得到后代    
%     end                                                     % find( ,1) 返回第一个 符合前面条件的数的 下标    
        
%     %状态估计，重采样以后，每个粒子的权重都变成了1/N    
%     x_est = mean(x_P);    
%         
%     % Save data in arrays for later plotting    
%     x_out = [x_out x];    
%     z_out = [z_out z];    
%     x_est_out = [x_est_out x_est];   
%     
    
% function outIndex = multinomialR(1:N;q)
% %获取数据长度
% outIndex=zeros(size(N));
% Col=length(q);
% N_babies= zeros(1,Col);
% 
% %计算粒子权重累计函数cdf 
% cdf= cumsum(q);
%  %产生[0,1]均匀分布的随机数
% u=rand(1,Col);
% 
% %求u^(j^-1)次方 
% uu=u.^(1./(Col:-1:1));
%  %如果A是一个向量，cumprod(A)将返回一个包含A各元素积累连乘的结果的向量
%  %元素个数与原向量相同
% ArrayTemp=cumprod(uu);
%  %fliplr(X)使矩阵X沿垂直轴左右翻转
% u = fliplr(ArrayTemp);
% j=1;
% for i=1:Col
%     %此处跟随机采样相似
%     while (u(i)>cdf(j))
%         j=j+1;
%     end
%     N_babies(j)=N_babies(j)+1;
% end;
% index=1;
% for i=1:Col
%     if (N_babies(i)>0)
%         for j=index:index+N_babies(i)-1
%             outIndex(j) = i;
%         end;
%     end;
%     index= index+N_babies(i);
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function outIndex = residualR(inIndex,q)
% outIndex=zeros(size(inIndex));
% N= length(q);
% N_babies= zeros(1,N);
% q_res = N.*q;
% N_babies = fix(q_res);
% N_res=N-sum(N_babies);
% if (N_res~=0)
%     q_res=(q_res-N_babies)/N_res;
%     cumDist= cumsum(q_res);
%     u = fliplr(cumprod(rand(1,N_res).^(1./(N_res:-1:1))));
%     j=1;
%     for i=1:N_res
%         while (u(1,i)>cumDist(1,j))
%             j=j+1;
%         end
%         N_babies(1,j)=N_babies(1,j)+1;
%     end;
% end;
% index=1;
% for i=1:N
%     if (N_babies(1,i)>0)
%         for j=index:index+N_babies(1,i)-1
%             outIndex(j) = i;
%         end;
%     end;
%     index= index+N_babies(1,i);
% end

 % 系统重采样
function outIndex = systematicR(inIndex,q); 
if nargin < 0, error('Not enough input arguments.'); end

q=q';
[arb,N] = size(q);   

 

N_children=zeros(1,N);
label=zeros(1,N);
label=1:1:N;

s=1/N;
auxw=0;
auxl=0;
li=0;   
T=s*rand(1);
j=1;
Q=0;
i=0;

 
u=rand(1,N);
while (T<1)
   if (Q>T)
      T=T+s;
      N_children(1,li)=N_children(1,li)+1;
   else
 
      i=fix((N-j+1)*u(1,j))+j;
 
      auxw=q(1,i);
      li=label(1,i);
 
      Q=Q+auxw;
 
      q(1,i)=q(1,j);
      label(1,i)=label(1,j);
   
      j=j+1;
   end
end
 
index=1;
for i=1:N
  if (N_children(1,i)>0)
    for j=index:index+N_children(1,i)-1
      outIndex(j) = inIndex(i);
    end;
  end;   
  index= index+N_children(1,i);   
end
