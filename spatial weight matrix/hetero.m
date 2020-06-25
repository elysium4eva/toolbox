function v=hetero(x,B,N)
gg=tril(ones(N),-1);
x1=x(1:N*(N-1)/2);
gg(gg==1)=x;
g=gg+gg';
sigma=diag(x(N*(N-1)/2+1:N*(N+1)/2));  
v=g*sigma^(-1)*g-B^(-1);