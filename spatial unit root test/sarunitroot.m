function [Y,d1,d2] = sarunitroot(N,m,n,T)
% Generate N panels with first m unit root panels for T periods
% m<=N
% using n heads and tails spatial dependence weight matrix
% 2n<N
% Y is the result of N by T matrix
% C is the cointegration matrix
W=head(N,n);
rho=rand(N,1);
T0=50;
Y(:,1)=zeros(N,1);
% alpha=normrnd(1,1,[1,N])';
% alpha=rand;
d2=eig(W);
rhow=(eye(N)-diag(rho)*W)^(-1);
lambda=zeros(N,1);
resid=zeros(N,1);

for i=1:N
lambda(i)=(1-rho(i))*rand;
resid(i)=1-lambda(i)-rho(i);
end
Resid=resid';
gamma=[Resid(1:m)'
    Resid(m+1:N)'*0.9];


for t=2:T+T0
   Y(:,t)=rhow*(diag(lambda)*Y(:,t-1)+diag(gamma)*W*Y(:,t-1)+normrnd(0,0.1,[N,1]));
end
Y=Y(:,T0+1:T+T0);
C=eye(N)-rhow*(diag(lambda)+diag(gamma)*W);
d1=eig(C);

%% 2SLS estimation of Y
WY=W*Y;
WYLL=W*Y(:,1:T-2);
WWYL=W*W*Y(:,2:T-1);
WWYLL=W*W*Y(:,1:T-2);
WWWYL=W*W*W*Y(:,2:T-1);
WWWWYL=W*W*W*W*Y(:,2:T-1);
beta=zeros(3,N);

for i=1:N
%% define the regressors
x=[WY(i,3:T)
   Y(i,2:T-1)
   WY(i,2:T-1)]';
y=Y(i,3:T)';
%% define the feasible instruments H
H=[Y(i,2:T-1)',WY(i,2:T-1)',WWYL(i,:)',WYLL(i,:)'];
%% define the projection matrix PH
PH=H*(H'*H)^(-1)*H';
beta(:,i)=(x'*PH*x)^(-1)*x'*PH*y;
end
rhowhat=(eye(N)-diag(beta(1,:))*W)^(-1);
Chat=eye(N)-rhowhat*(diag(beta(2,:))+diag(beta(3,:))*W);
d3=eig(Chat);
