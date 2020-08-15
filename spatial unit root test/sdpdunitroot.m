function [results] = sdpdunitroot(N,n,T,rho,info)
% Generate N panels with first m unit root panels for T periods
% m<=N
% using n heads and tails spatial dependence weight matrix
% 2n<N
% Y is the result of N by T matrix
% C is the cointegration matrix
W=head(N,n);
T0=50;
Y(:,1)=zeros(N,1);
rhow=(eye(N)-rho*W)^(-1);
lambda=zeros(N,1);
gamma=zeros(N,1);

for i=1:N
lambda(i)=(1-rho)*rand;
gamma(i)=(1-lambda(i)-rho)*rand;
end

if info.model==0
for t=2:T+T0
   Y(:,t)=rhow*(diag(lambda)*Y(:,t-1)+diag(gamma)*W*Y(:,t-1)+normrnd(0,0.1,[N,1]));
end

elseif info.model==1
   mu=rand(N,1)*0.2;
   for t=2:T+T0
   Y(:,t)=rhow*(diag(lambda)*Y(:,t-1)+diag(gamma)*W*Y(:,t-1)+mu+normrnd(0,0.1,[N,1]));
   end 
else
    mu=rand(N,1)*0.2;
    delta=rand(N,1)*0.02;
       for t=2:T+T0
   Y(:,t)=rhow*(diag(lambda)*Y(:,t-1)+diag(gamma)*W*Y(:,t-1)+mu+delta*t+normrnd(0,1,[N,1]));
       end 
end
Y=Y(:,T0+1:T+T0);
C=eye(N)-rhow*(diag(lambda)+diag(gamma)*W);
d1=eig(C);

%% 2SLS estimation of Y
results=sdpd(Y,W,info);



