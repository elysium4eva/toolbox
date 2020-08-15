% Produces column (6) of Tables 4 in the SLX paper
%
A=xlsread('x:\lotus\cigarette+2var.xls');
B=xlsread('x:\lotus\cigar_states.xls');
%W1=wk1read('x:\lotus\Spat-Sym-US.wk1');
W1=xlsread('x:\lotus\Spat-Sym-US.xls');
T=30; % number of time periods
N=46; % number of regions
W=normw(W1); % function of LeSage
y=A(1:end,3); % column number in the data matrix that corresponds to the dependent variable
x=A(1:end,[4,6]); % column numbers in the data matrix that correspond to the independent variables (4 and 6)
K=4;
pop=A(1:end,8)/1000;
tax=A(1:end,9);
comp=A(1:end,10)/1000;
for i=1:N
    for j=1:N
        if (i==j) 
            W(i,j)=0;
        else
            W(i,j)=1/sqrt((B(i,1)-B(j,1))^2+(B(i,2)-B(j,2))^2); %Pythagoras
        end
    end
end
lambda=eig(W);
W=W/max(lambda); 
for t=1:T
    t1=(t-1)*N+1;t2=t*N;
    wpop(t1:t2,1)=W*pop(t1:t2,1); 
    wtax(t1:t2,1)=W*tax(t1:t2,1);
    wcomp(t1:t2,1)=W*comp(t1:t2,1);
end

% Endogeneity test SLX model with W=Parameterized inverse distance
% First exogenous estimation

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
for t=1:T
    t1=(t-1)*N+1;t2=t*N;
    wx(t1:t2,:)=Walpha*x(t1:t2,:);
end
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,[x wx],N,T,model);
results=ols(ywith,xwith);
vnames=strvcat('logcit','logp','logy','W*logp','W*logy');
prt_reg(results,vnames);
loglikelihood = -(N*T/2)*log(2*pi) - (N*T/2)*log(results.si2) - (1/results.si2)*results.resid'*results.resid
parm=[results.beta;alpha];
xpx = hessian('f2_gravity_alpha',parm,y,x,W,T,model); 
[ndummy nvar]=size(xwith);
xpx(1:nvar,1:nvar)=xwith'*xwith; 
varcov=results.sige*invpd(xpx); 
tvalues=parm./sqrt(abs(diag(varcov)));
[parm tvalues]

fprintf(1,'reprinted, W = 1/d^alpha, alpha estimated, also normalized by largest eigenvalue\n'); %normalization ensures, I believe, more sensible magnitudes of the coefficient estimates wx1 and wx2
lambda=eig(Walpha);
WT=Walpha/max(lambda);
for t=1:T
    t1=(t-1)*N+1;t2=t*N;
    wx(t1:t2,:)=WT*x(t1:t2,:);
end
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,[x wx],N,T,model);
results=ols(ywith,xwith);
prt_reg(results,vnames);
loglikelihood = -(N*T/2)*log(2*pi) - (N*T/2)*log(results.si2) - (1/results.si2)*results.resid'*results.resid
parm=[results.beta;alpha];
fprintf(1,'parameter estimates and t-values including alpha\n');
[parm tvalues]

intercept=mean(y)-mean([x wx])*results.beta; 
sfe=meanny-meannx*results.beta-kron(en,intercept); %state fixed effects  vector 46*1
tfe=meanty-meantx*results.beta-kron(et,intercept);
yme = y - mean(y); % vector 1380*1
ent=ones(N*T,1);
error=y-kron(tfe,en)-kron(et,sfe)-[x wx]*results.beta-kron(ent,intercept);
rsqr1 = error'*error;
rsqr2 = yme'*yme;
FE_rsqr2 = 1.0 - rsqr1/rsqr2 % r-squared including fixed effects
sige=results.sige*((N*T-K)/(N*T));
loglik=-N*T/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid


% Endogeneity test preparation
for t=1:T
    wpop(t1:t2,1)=WT*pop(t1:t2,1); 
    wtax(t1:t2,1)=WT*tax(t1:t2,1);
    wcomp(t1:t2,1)=WT*comp(t1:t2,1);
end
vari=[pop wtax wcomp];

%for k=1:3
%
%k
Ninstr=2;
%xinstr=[pop wpop tax wtax comp wcomp x(:,2) wx(:,2)]; % population + W, tax + W, comp + w, income + W, No lagged consumption
%xinstr=[wpop tax comp x(:,2)]; % population + W, tax + W, comp + w, income + W, No lagged consumption
%xinstr=[vari(:,k) tax wpop wcomp x(:,2) wx(:,2)]; % population + W, tax + W, comp + w, income + W, No lagged consumption
xinstr=[tax wpop wcomp x(:,2) wx(:,2)]; % population + W, tax + W, comp + w, income + W, No lagged consumption
%Nx=8;
%Nx=4;
Nx=5;
Noverinstr=Nx-Ninstr;

% First potential endogenous regressor
yinstr1=x(:,1);
[yinstrwith1,xinstrwith,meanny,meannx,meanty,meantx]=demean(yinstr1,xinstr,N,T,3);
corrcoef(xinstrwith)
results=ols(yinstrwith1,xinstrwith);
% Correction for degrees of freedom
results.sige=results.sige*(N*T-Nx)/(N*T-Nx-N-(T-1)); 
tmp = (results.sige)*(diag(results.xpxi));
results.tstat = results.beta./(sqrt(tmp));
% end correction
vnames=strvcat('logp','tax','wpop','wcomp','logy','Wlogy');
%vnames=strvcat('logp','wpop','tax','comp','logy');
%vnames=strvcat('logp','pop','wpop','tax','wtax','comp','wcomp','logy','Wlogy');
prt_reg(results,vnames);
errorend1=results.resid;
% F-test on significance instruments
btemp=results.beta(1:Nx);
varcov=results.sige*results.xpxi(1:Nx,1:Nx);
Waldtest=btemp'*inv(varcov)*btemp;
prob_Waldtest=1-chis_cdf(Waldtest,Nx); % probability greater than 0.05 points to insignificance
[Waldtest prob_Waldtest]
sige=results.sige*(N*T-Nx)/(N*T-Nx-N-(T-1));
Ftest=Waldtest/Nx;
prob_Ftest=1-fdis_cdf(Ftest,Nx,N*T-Nx-N-(T-1)); % probability greater than 0.05 points to insignificance
[Ftest prob_Ftest]

% Second potential endogenous regressor
yinstr2=wx(:,1);
[yinstrwith2,xinstrwith,meanny,meannx,meanty,meantx]=demean(yinstr2,xinstr,N,T,3);
results=ols(yinstrwith2,xinstrwith);
% Correction for degrees of freedom
results.sige=results.sige*(N*T-Nx)/(N*T-Nx-N-(T-1)); 
tmp = (results.sige)*(diag(results.xpxi));
results.tstat = results.beta./(sqrt(tmp));
% end correction
vnames=strvcat('W*logp','tax','wpop','wcomp','logy','Wlogy');
%vnames=strvcat('W*logp','wpop','tax','comp','logy');
%vnames=strvcat('W*logp','pop','wpop','tax','wtax','comp','wcomp','logy','Wlogy');
prt_reg(results,vnames);
errorend2=results.resid;
% F-test on significance instruments
btemp=results.beta(1:Nx);
varcov=results.sige*results.xpxi(1:Nx,1:Nx);
Waldtest=btemp'*inv(varcov)*btemp;
prob_Waldtest=1-chis_cdf(Waldtest,Nx); % probability greater than 0.05 points to insignificance
[Waldtest prob_Waldtest]
sige=results.sige*(N*T-Nx)/(N*T-Nx-N-(T-1));
Ftest=Waldtest/Nx;
prob_Ftest=1-fdis_cdf(Ftest,Nx,N*T-Nx-N-(T-1)); % probability greater than 0.05 points to insignificance
[Ftest prob_Ftest]

% OLS model with residual IV regression
fprintf(1,'OLS with residuals IV regression\n');
results=ols(ywith,[errorend1 errorend2 xwith]); % 
vnames=strvcat('logcit','resid1 IV','resid2 IV','logp','logy','W*logp','W*logy');
% Correction for degrees of freedom
results.sige=results.sige*(N*T-K-1)/(N*T-K-1-N-(T-1)); 
tmp = (results.sige)*(diag(results.xpxi));
results.tstat = results.beta./(sqrt(tmp));
% end correction
prt_reg(results,vnames);

% 2SLS regression
fprintf(1,'2SLS model\n');
yinstrwith=[yinstrwith1 yinstrwith2];
results=tsls(ywith,yinstrwith,xwith(:,[2,4]),xinstrwith);
vnames=strvcat('logcit','logp','W*logp','logy','W*logy');
prt_reg(results,vnames);
resiv=results.resid;
% Validity check IV
fprintf(1,'Validity IV\n');
results=ols(resiv,xinstrwith); % 
vnames=strvcat('resiv','tax','wpop','wcomp','logy','Wlogy');
%vnames=strvcat('resiv','wpop','tax','comp','logy');
%vnames=strvcat('resiv','pop','wpop','tax','wtax','comp','wcomp','logy','Wlogy');
% Correction for degrees of freedom
results.sige=results.sige*(N*T-K)/(N*T-K-N-(T-1)); 
tmp = (results.sige)*(diag(results.xpxi));
results.tstat = results.beta./(sqrt(tmp));
% end correction
prt_reg(results,vnames);
Vtest=N*T*results.rsqr;
prob_Vtest=1-chis_cdf(Vtest,Noverinstr); % probability greater than 0.05 points to insignificance
fprintf(1,'Sargan\n');
[Vtest prob_Vtest]

%end


% 2SLS estimation problem: if W changes the endogenous variable WX changes,
% this complicates the estimation
% 
alphaold=1;
iter=0;converge=1.0;criteria=0.0001;itermax=25;
options.Display='off';
options.MaxFunEvals=1000;
options.MaxIter=1000;
options.TolX=0.000001;
options.TolFun=0.000001;
model=3;
yinstr=[x(:,1) wx(:,1)]; % endogenous explanatary variable(s)
xi=x(:,2); % exogenous explanatory variable(s)
wxi=x(:,2); % SLX exogenous explanatory variable(s), to be multiplied by W
xinstr=[tax x(:,2)]; % first set of instruments
wxinstr=[pop comp x(:,2)]; % second set of instruments, to be multiplied by W
[alpha fval exitflag]=fminsearch('f_gravity_endogenous_alpha',alphaold,options,y,yinstr,xi,wxi,xinstr,wxinstr,W,T,model) 
Walpha=W.^alpha;
for t=1:T
    t1=(t-1)*N+1;t2=t*N;
    Wxi(t1:t2,:)=Walpha*wxi(t1:t2,:);
    Wxinstr(t1:t2,:)=Walpha*wxinstr(t1:t2,:);
end
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,[xi Wxi],N,T,model);
[ywith,yinstrwith,meanny,meannx,meanty,meantx]=demean(y,yinstr,N,T,model);
[ywith,xinstrwith,meanny,meannx,meanty,meantx]=demean(y,[xinstr Wxinstr],N,T,model);
results=tsls(ywith,yinstrwith,xwith,xinstrwith);
vnames=strvcat('logcit','logp','logy','W*logp','W*logy');
prt_reg(results,vnames);

parm=[results.beta;alpha];
xpx = hessian('f2_gravity_endogenous_alpha',parm,y,yinstr,xi,wxi,xinstr,wxinstr,W,T,model);
[ndummy nvar]=size([yinstrwith xwith]);
xpx(1:nvar,1:nvar)=results.xpx; 
varcov=results.sige*invpd(xpx); 
tvalues=parm./sqrt(abs(diag(varcov)));
lambda=eig(Walpha);
WT=Walpha/max(lambda);
for t=1:T
    t1=(t-1)*N+1;t2=t*N;
    Wxi(t1:t2,:)=WT*wxi(t1:t2,:);
    Wxinstr(t1:t2,:)=WT*wxinstr(t1:t2,:);
end
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,[xi Wxi],N,T,model);
[ywith,yinstrwith,meanny,meannx,meanty,meantx]=demean(y,yinstr,N,T,model);
[ywith,xinstrwith,meanny,meannx,meanty,meantx]=demean(y,[xinstr Wxinstr],N,T,model);
results=tsls(ywith,yinstrwith,xwith,xinstrwith);
vnames=strvcat('logcit','logp','W*logp','logy','W*logy');
prt_reg(results,vnames);

[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,[x wx],N,T,model);
intercept=mean(y)-mean([x wx])*results.beta; 
sfe=meanny-meannx*results.beta-kron(en,intercept); %state fixed effects  vector 46*1
tfe=meanty-meantx*results.beta-kron(et,intercept);
yme = y - mean(y); % vector 1380*1
ent=ones(N*T,1);
error=y-kron(tfe,en)-kron(et,sfe)-[x wx]*results.beta-kron(ent,intercept);
rsqr1 = error'*error;
rsqr2 = yme'*yme;
FE_rsqr2 = 1.0 - rsqr1/rsqr2 % r-squared including fixed effects
parm=[results.beta;alpha];
[parm tvalues]