%  程序说明： 粒子滤波子程序
function [Xout,Tpf]=PF(Z,NodePostion,canshu)
  
M=canshu.M;
Q=canshu.Q;
R=canshu.R;
F=canshu.F;
T=canshu.T;
state0=canshu.state0;
Node_number=canshu.Node_number;
 
N=300; 
r=zeros(1,N);
zPred=zeros(1,N);
Weight=zeros(1,N);
xparticlePred=zeros(4,N);
Xout=zeros(4,M);
Xout(:,1)=state0;
Tpf=zeros(1,M);
for i=1:Node_number
    xparticle{i}=zeros(4,N);
    for j=1:N    
        xparticle{i}(:,j)=state0;
    end
    Xpf{i}=zeros(4,N);
    Xpf{i}(:,1)=state0;
end
 
for t=2:M
    tic;
    XX=0;
    for i=1:Node_number
        x0=NodePostion(1,i);
        y0=NodePostion(2,i);
 
        for k=1:N
            xparticlePred(:,k)=feval('sfun',xparticle{i}(:,k),T,F)+5*sqrtm(Q)*randn(4,1);
        end
 
        for k=1:N
            zPred(1,k)=feval('hfun',xparticlePred(:,k),x0,y0);
        end
        for k=1:N
            r(1,k)=Z(1,t)-zPred(1,k);
        end

        XzPred(1,:)=PSO(zPred(1,:),r(1,:));
        
        for k=1:N
            z1=Z(i,t)-XzPred(1,k);
            Weight(1,k)=inv(sqrt(2*pi*det(R)))*exp(-.5*(z1)'*inv(R)*(z1))+ 1e-99; 
        end
        
        Weight(1,:)=Weight(1,:)./sum(Weight(1,:));
        outIndex = systematicR(1:N,Weight(1,:)');        
        xparticle{i}= xparticlePred(:,outIndex);  
        target=[mean(xparticle{i}(1,:)),mean(xparticle{i}(2,:)),...
            mean(xparticle{i}(3,:)),mean(xparticle{i}(4,:))]';
        Xpf{i}(:,t)=target;
      
        XX=XX+Xpf{i}(:,t);
    end
    Xout(:,t)=XX/Node_number;
%   xparticle{i}= xparticlePred(:,Weight(1,:));
    Tpf(1,t)=toc;
end

% %  程序说明： 粒子滤波子程序
% function [Xout,Tpf]=PF(Z,NodePostion,canshu)
% M=canshu.M;
% Q=canshu.Q;
% R=canshu.R;
% F=canshu.F;
% T=canshu.T;
% state0=canshu.state0;
% Node_number=canshu.Node_number;
%  
% N=100;       
% zPred=zeros(1,N);
% Weight=zeros(1,N);
% xparticlePred=zeros(4,N);
% Xout=zeros(4,M);
% Xout(:,1)=state0;
% Tpf=zeros(1,M);
% for i=1:Node_number
%     xparticle{i}=zeros(4,N);
%     for j=1:N    
%         xparticle{i}(:,j)=state0;
%     end
%     Xpf{i}=zeros(4,N);
%     Xpf{i}(:,1)=state0;
% end
%  
% for t=2:M
%     tic;
%     XX=0;
%     for i=1:Node_number
%         x0=NodePostion(1,i);
%         y0=NodePostion(2,i);
%  
%         for k=1:N
%             xparticlePred(:,k)=feval('sfun',xparticle{i}(:,k),T,F)+5*sqrtm(Q)*randn(4,1);
% %             xpartminus(i) = 0.5 * xpart(i) + 25 * xpart(i) / (1 + xpart(i)^2) ...
% %          + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
% %             xparticlePred(:,k)=0.5*xparticle{i}(:,k)+25*xparticle{i}(:,k)/(xparticle{i}(:,k)*(xparticle{i}(:,k))')+sqrt(Q)*randn(4,1);
%         end
% %  
%         for k=1:N
%             zPred(1,k)=feval('hfun',xparticlePred(:,k),x0,y0);
%             z1=Z(i,t)-zPred(1,k);
%             Weight(1,k)=inv(sqrt(2*pi*det(R)))*exp(-.5*(z1)'*inv(R)*(z1))+ 1e-99; 
% %             q(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R); 
% %             Weight(1,k)=(1 / sqrt(R) / sqrt(2*pi)) * exp(-z1^2 / 2 / R); 
%         end
%  
%         Weight(1,:)=Weight(1,:)./sum(Weight(1,:));
%         outIndex = randomR(1:N,Weight(1,:)');        
%         xparticle{i}= xparticlePred(:,outIndex);  
%         target=[mean(xparticle{i}(1,:)),mean(xparticle{i}(2,:)),...
%             mean(xparticle{i}(3,:)),mean(xparticle{i}(4,:))];
%         Xpf{i}(:,t)=target;
%       
%         XX=XX+Xpf{i}(:,t);
%     end
%     Xout(:,t)=XX/Node_number;
%     Tpf(1,t)=toc;
%     save 目标预测.txt xparticlePred -ascii;
%     save 真实量测.txt zPred -ascii;
% end
