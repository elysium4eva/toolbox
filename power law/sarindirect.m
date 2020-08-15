function [results] = sarindirect(rho,beta,w,N,T)
% calculate degree distribution of sar model
% rho is the estimated SAR coefficient
% w is the prespecified weight matrix
wt=zeros(N,N);
invwt=zeros(N,N);
results.d1(1,:)=zeros(1,N);
results.d2(1,:)=zeros(1,N);

for t=1:T
    wt=w(((t-1)*N+1:t*N),((t-1)*N+1:t*N));
    invwt=(eye(N)-rho*wt)^(-1);
    results.d1(t,:)=sum(invwt-diag(diag(invwt)))*beta; 
    results.d2(t,:)=diag(invwt)*beta;
end

