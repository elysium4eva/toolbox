% Endogeneity test preparation
Ninstr=1;
%xinstr=[pop wpop tax wtax comp wcomp x(:,2) wx(:,2)];
%vari=[pop wpop wtax comp wcomp];
%for k=1:5
%k 
xinstr=[lnapp wlntfp lnhr wlnhr wwlns];
Nx=5;
Noverinstr=Nx-Ninstr;
% Potential endogenous regressor
yinstr1=wlns;
[yinstrwith1,xinstrwith,meanny,meannx,meanty,meantx]=demean(yinstr1,xinstr,N,T,1);
results=ols(yinstrwith1,xinstrwith);
% Correction for degrees of freedom
results.sige=results.sige*(N*T-Nx)/(N*T-Nx-N-(T-1)); 
tmp = (results.sige)*(diag(results.xpxi));
results.tstat = results.beta./(sqrt(tmp));
% end correction
% vnames=strvcat('logp','pop','wpop','tax','wtax','comp','wcomp','logy','Wlogy');
vnames=strvcat('wlns','lnapp','wlntfp','lnhr','wlnhr','wwlns');
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

%% fixed effect demean of dependent and independent variables
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x2,N,T,1);
K=3;

% OLS model with residual IV regression
fprintf(1,'OLS with residuals IV regression\n');
results=ols(ywith,[errorend1 xwith]); % 
vnames=strvcat('logtfp','resid1 IV','lnhr','lns');
% Correction for degrees of freedom
results.sige=results.sige*(N*T-K-1)/(N*T-K-1-N-(T-1)); 
tmp = (results.sige)*(diag(results.xpxi));
results.tstat = results.beta./(sqrt(tmp));
% end correction
prt_reg(results,vnames);

% 2SLS regression
fprintf(1,'2SLS model\n');
yinstrwith=[yinstrwith1];
results=tsls(ywith,yinstrwith,xwith,xinstrwith);
vnames=strvcat('logtfp','wlns','lnhr','lns');
prt_reg(results,vnames);
resiv=results.resid;
% Validity check IV
fprintf(1,'Validity IV\n');
results=ols(resiv,xinstrwith); % 
vnames=strvcat('resiv','lnapp','wlntfp','lnhr','wlnhr','wwlns');
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