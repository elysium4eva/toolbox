function [d] = sarindirect(rho,beta,w,N,T)
% calculate degree distribution of sar model
% rho is the estimated SAR coefficient
% w is the prespecified weight matrix
wt=zeros(N,N);
invwt=zeros(N,N);
d(1,:)=zeros(1,N);

for t=1:T
    wt=w(((t-1)*N+1:t*N),((t-1)*N+1:t*N));
    invwt=(eye(N)-rho*wt)^(-1)*beta;
    d(t,:)=sum(invwt); 
end

