function llike=f_gravity_alpha(alpha,y,x,W,T,model)
Walpha=W.^alpha;
[N N1]=size(W);
for t=1:T
    t1=(t-1)*N+1;t2=t*N;
    wx(t1:t2,:)=Walpha*x(t1:t2,:);
end
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,[x wx],N,T,model);
results=ols(ywith,xwith);
res=results.resid;
llike=res'*res;