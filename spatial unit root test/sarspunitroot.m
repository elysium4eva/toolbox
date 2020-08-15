function [Y,d1,d2] = sarspunitroot(N,n,T)
% Generate N panels with first m unit root panels for T periods
% m<=N
% using n heads and tails spatial dependence weight matrix
% 2n<N
% Y is the result of N by T matrix
% C is the cointegration matrix
W=head(N,n);
rho=rand(N,1);
Y0(:,1)=zeros(N,1);
% alpha=normrnd(1,1,[1,N])';
% alpha=rand;
d2=eig(W);
rhow=(eye(N)-diag(rho)*W)^(-1);

lambda=ones(N,1)-rho;
gamma=0.9*rho;


for t=2:T
   for i=1:N
     %  Y0(i,t)=alpha+gamma(i)*Y0(i,t-1)+lambda(i)*W(i,:)*Y0(:,t-1)+normrnd(0,1);
     %  Y0(i,t)=alpha(i)+gamma(i)*Y0(i,t-1)+lambda(i)*W(i,:)*Y0(:,t-1)+normrnd(0,1);
     Y0(i,t)=gamma(i)*Y0(i,t-1)+lambda(i)*W(i,:)*Y0(:,t-1)+normrnd(0,1);
   end
end
Y=rhow*Y0;
C=eye(N)-diag(gamma)-diag(lambda)*W;
d1=eig(C);

end