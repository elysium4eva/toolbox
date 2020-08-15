function llike = f_sec1(parm,uhat,W,N,T)
rsqr=0;
for t=1:T
    u=uhat(:,t);
    epsilon=(eye(N)-parm*W)*u;
    rsqr=rsqr+epsilon'*epsilon;
end
e=eig(W);
lndet=0;
for i=1:N
    lndet=lndet+log(1-e(i)*parm);    
end
llike=0.5*rsqr-(T)*lndet;


