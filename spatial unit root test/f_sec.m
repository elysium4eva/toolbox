function llike = f_sec(parm,y,W,N,T,T0)
rsqr=0;
for t=1:T+T0
    u=y(:,t+1)-y(:,t);
    epsilon=(eye(N)-parm*W)*u;
    rsqr=rsqr+epsilon'*epsilon;
end
e=eig(W);
lndet=0;
for i=1:N
    lndet=lndet+log(1-e(i)*parm);    
end
llike=0.5*rsqr-(T+T0)*lndet;

