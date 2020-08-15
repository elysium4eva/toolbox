function llike = f_sdpd(parm,y,W,info)
rsqr=0;
Wy=W*y;
[N,T]=size(y);
trend=1:T-1;
if info.model==0
for i=1:N
    x=[y(i,1:T-1)' Wy(i,1:T-1)'];
    Mx=eye(T-1)-x*(x'*x)^(-1)*x';
    sy=y(i,2:T)-parm*Wy(i,2:T);
    rsqr=rsqr+log(sy*Mx*sy'/T);
end
elseif info.model==1
    for i=1:N
    x=[ones(T-1,1) y(i,1:T-1)' Wy(i,1:T-1)'];
    Mx=eye(T-1)-x*(x'*x)^(-1)*x';
    sy=y(i,2:T)-parm*Wy(i,2:T);
    rsqr=rsqr+log(sy*Mx*sy'/T);
    end
else
    for i=1:N
    x=[ones(T-1,1) trend' y(i,1:T-1)' Wy(i,1:T-1)'];
    Mx=eye(T-1)-x*(x'*x)^(-1)*x';
    sy=y(i,2:T)-parm*Wy(i,2:T);
    rsqr=rsqr+log(sy*Mx*sy'/T);
    end
    
end
e=eig(W);
lndet=0;
for i=1:N
    lndet=lndet+log(1-e(i)*parm);    
end
llike=0.5*T*rsqr-T*lndet;


