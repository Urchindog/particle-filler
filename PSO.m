function outIndex = PSO(inIndex,r)
%% ��ʼ����Ⱥ
% f=inv(sqrt(2*pi*det(R)))*exp(-.5*(z1)'*inv(R)*(z1))+ 1e-99; 
% f= @(x)x .* sin(x) .* cos(2 * x) - 2 * x .* sin(3 * x); % �������ʽ
% figure(1);ezplot(f,[0,0.01,20]);
f=r;
N = 300;                         % ��ʼ��Ⱥ����
d = 1;                          % �ռ�ά��
ger = 1000;                      % ����������     
limit = [0, 100];               % ����λ�ò�������
vlimit = [-3, 3];               % �����ٶ�����
w = 0.8;                        % ����Ȩ��
c1 = 1.5;                       % ����ѧϰ����
c2 = 1.5;                       % Ⱥ��ѧϰ���� 
x=zeros(1,N);
% for i = 1:d
%     x = limit(i, 1) + (limit(i, 2) - limit(i, 1)) * rand(N, d);%��ʼ��Ⱥ��λ��
% end
x=inIndex;
v = rand(d,N);                  % ��ʼ��Ⱥ���ٶ�
xm = zeros(1, N);                 % ÿ���������ʷ���λ��
ym = zeros(1, d);                % ��Ⱥ����ʷ���λ��
fxm = zeros(1,N);               % ÿ���������ʷ�����Ӧ��
fym = -inf;                      % ��Ⱥ��ʷ�����Ӧ��

%% Ⱥ�����
iter = 1;
record = zeros(ger, 1);          % ��¼��
while iter <= ger
     fx = r ;                 % ���嵱ǰ��Ӧ��   
     for g = 1:N      
        if fxm(g) < fx(g)
            fxm(g) = fx(g);      % ���¸�����ʷ�����Ӧ��
            xm(1,g) = x(1,g);    % ���¸�����ʷ���λ��
        end 
     end
if fym < min(fxm)
        fym = min(fxm);   % ����Ⱥ����ʷ�����Ӧ��
        ym = min(xm);         % ����Ⱥ����ʷ���λ��
 end
    v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, 1, N) - x);% �ٶȸ���
    % �߽��ٶȴ���
    v(v > vlimit(2)) = vlimit(2);
    v(v < vlimit(1)) = vlimit(1);
    x = x + v;% λ�ø���
    % �߽�λ�ô���
    x(x > limit(2)) = limit(2);
    x(x < limit(1)) = limit(1);
    record(iter) = fym;%���ֵ��¼
%     x0 = 0 : 0.01 : 20;
%     plot(x0, f(x0), 'b-', x, f(x), 'ro');title('״̬λ�ñ仯')
%     pause(0.1)
    iter = iter+1;
end
outIndex=x;



% figure(3);plot(record);title('��������')
% x0 = 0 : 0.01 : 20;
% figure(4);plot(x0, f(x0), 'b-', x, f(x), 'ro');title('����״̬λ��')
% disp(['���ֵ��',num2str(fym)]);
% disp(['����ȡֵ��',num2str(ym)]);
