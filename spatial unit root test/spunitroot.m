function [Y,d1,d2,W,rho] = spunitroot(N,m,n,T)
% Generate N panels with first m unit root panels for T periods
% m<=N
% using n heads and tails spatial dependence weight matrix
% 2n<N
% Y is the result of N by T matrix
% C is the cointegration matrix
W=head(N,n);
rho=rand;

Y(:,1)=zeros(N,1);
% alpha=normrnd(1,1,[1,N])';
% alpha=rand;
d2=eig(W);
rhow=(eye(N)-rho*W)^(-1);


Gamma=ones(N,1)*(1-rho);
gamma=[Gamma(1:m,1)
    Gamma(m+1:N,1)*0.9];


for t=2:T
   Y(:,t)=rhow*(diag(gamma)*W*Y(:,t-1)+normrnd(0,1,[N,1]));
end
C=eye(N)-rhow*(diag(gamma)*W);
d1=eig(C);

end