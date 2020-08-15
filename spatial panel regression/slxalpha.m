function [results] = slxalpha(y,x,w,N,T,vnames)
%OLSALPHA estimate the slx model with parameterized alpha, alpha is the
%exponent of inverse distance
%Fixed effect models is estimated
%% vnames is variable name for  "y" "x" "wx" to display the resutls
lambda=eig(w);
W=w/max(lambda); 
[nob,k]=size(x);
K=2*k;

fprintf(1,'OLS model with exogenous interaction effects, W = Parameterized inverse distance\n');

% Estimating one alpha
model=3;
en=ones(N,1);
et=ones(T,1);
alphaold=1;
iter=0;converge=1.0;criteria=0.0001;itermax=25;
options.Display='off';
options.MaxFunEvals=1000;
options.MaxIter=1000;
options.TolX=0.000001;
options.TolFun=0.000001;

[alpha fval exitflag]=fminsearch('f_gravity_alpha',alphaold,options,y,x,W,T,model); 

Walpha=W.^alpha;
lambda=eig(Walpha);
WT=Walpha/max(lambda);
for t=1:T
    t1=(t-1)*N+1;t2=t*N;
    wx(t1:t2,:)=WT*x(t1:t2,:);
end
fprintf(1,'reprinted, W = 1/d^alpha, alpha estimated, also normalized by largest eigenvalue\n'); %normalization ensures, I believe, more sensible magnitudes of the coefficient estimates wx1 and wx2
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,[x wx],N,T,model);
results=ols(ywith,xwith);
prt_reg(results,vnames);
results.loglikelihood = -(N*T/2)*log(2*pi) - (N*T/2)*log(results.si2) - (1/results.si2)*results.resid'*results.resid;
results.parm=[results.beta;alpha];
xpx = hessian('f2_gravity_alpha',results.parm,y,x,W,T,model); 
[ndummy nvar]=size(xwith);
xpx(1:nvar,1:nvar)=xwith'*xwith; 
varcov=results.sige*invpd(xpx); 
tvalues=results.parm./sqrt(abs(diag(varcov)));
fprintf(1,'parameter estimates and t-values including alpha\n');
result.tstats=[results.parm tvalues];


results.intercept=mean(y)-mean([x wx])*results.beta; 
results.sfe=meanny-meannx*results.beta-kron(en,results.intercept); %region fixed effects  vector N by 1
results.tfe=meanty-meantx*results.beta-kron(et,results.intercept);
yme = y - mean(y); % vector N*T by 1
ent=ones(N*T,1);
error=y-kron(results.tfe,en)-kron(et,results.sfe)-[x wx]*results.beta-kron(ent,results.intercept);
rsqr1 = error'*error;
rsqr2 = yme'*yme;
results.FE_rsqr2 = 1.0 - rsqr1/rsqr2; % r-squared including fixed effects
sige=results.sige*((N*T-K)/(N*T));
results.FE_loglik=-N*T/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid;


end

