function llike=f2_gravity_alpha(parm,y,x,W,T,model)
beta=parm(1:end-1);
alpha=parm(end);
[N N1]=size(W);
Walpha=W.^alpha;
for t=1:T
    t1=(t-1)*N+1;t2=t*N;
    wx(t1:t2,:)=Walpha*x(t1:t2,:);
end
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,[x wx],N,T,model);
res=ywith-xwith*beta;
si2=1/(N*T)*res'*res;
%llike=0.5*res'*res;
llike=(N*T/2)*log(2*pi) + (N*T/2)*log(si2) + (1/si2)*res'*res;