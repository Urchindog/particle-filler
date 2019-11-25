% 程序说明： 观测方程函数
% 输入参数： x目标的状态，（x0,y0)是观测站的位置
% 输出参数： y是角度
function [y]=hfun(x,x0,y0)
 
if nargin < 3
    error('Not enough input arguments.'); 
end
[row,col]=size(x);
if row~=4|col~=1
    error('Input arguments error!');
end
xx=x(1)-x0;
yy=x(3)-y0;
y=atan2(yy,xx);
