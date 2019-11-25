function outIndex = PSO(inIndex,r)
%% 初始化种群
% f=inv(sqrt(2*pi*det(R)))*exp(-.5*(z1)'*inv(R)*(z1))+ 1e-99; 
% f= @(x)x .* sin(x) .* cos(2 * x) - 2 * x .* sin(3 * x); % 函数表达式
% figure(1);ezplot(f,[0,0.01,20]);
f=r;
N = 300;                         % 初始种群个数
d = 1;                          % 空间维数
ger = 1000;                      % 最大迭代次数     
limit = [0, 100];               % 设置位置参数限制
vlimit = [-3, 3];               % 设置速度限制
w = 0.8;                        % 惯性权重
c1 = 1.5;                       % 自我学习因子
c2 = 1.5;                       % 群体学习因子 
x=zeros(1,N);
% for i = 1:d
%     x = limit(i, 1) + (limit(i, 2) - limit(i, 1)) * rand(N, d);%初始种群的位置
% end
x=inIndex;
v = rand(d,N);                  % 初始种群的速度
xm = zeros(1, N);                 % 每个个体的历史最佳位置
ym = zeros(1, d);                % 种群的历史最佳位置
fxm = zeros(1,N);               % 每个个体的历史最佳适应度
fym = -inf;                      % 种群历史最佳适应度

%% 群体更新
iter = 1;
record = zeros(ger, 1);          % 记录器
while iter <= ger
     fx = r ;                 % 个体当前适应度   
     for g = 1:N      
        if fxm(g) < fx(g)
            fxm(g) = fx(g);      % 更新个体历史最佳适应度
            xm(1,g) = x(1,g);    % 更新个体历史最佳位置
        end 
     end
if fym < min(fxm)
        fym = min(fxm);   % 更新群体历史最佳适应度
        ym = min(xm);         % 更新群体历史最佳位置
 end
    v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, 1, N) - x);% 速度更新
    % 边界速度处理
    v(v > vlimit(2)) = vlimit(2);
    v(v < vlimit(1)) = vlimit(1);
    x = x + v;% 位置更新
    % 边界位置处理
    x(x > limit(2)) = limit(2);
    x(x < limit(1)) = limit(1);
    record(iter) = fym;%最大值记录
%     x0 = 0 : 0.01 : 20;
%     plot(x0, f(x0), 'b-', x, f(x), 'ro');title('状态位置变化')
%     pause(0.1)
    iter = iter+1;
end
outIndex=x;



% figure(3);plot(record);title('收敛过程')
% x0 = 0 : 0.01 : 20;
% figure(4);plot(x0, f(x0), 'b-', x, f(x), 'ro');title('最终状态位置')
% disp(['最大值：',num2str(fym)]);
% disp(['变量取值：',num2str(ym)]);
