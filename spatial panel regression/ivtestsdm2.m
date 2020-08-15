% Endogeneity test preparation
Ninstr=2;
%xinstr=[pop wpop tax wtax comp wcomp x(:,2) wx(:,2)];
%vari=[pop wpop wtax comp wcomp];
%for k=1:5
%k 
%xinstr=[lnapp lns wlns lnhr wlnhr];
xinstr=[wlns1 wwlns1 wlnhr1 wwlnhr1  lnapp];
[nobs,Nx]=size(xinstr);
Noverinstr=Nx-Ninstr;
% Potential endogenous regressor 1
yinstr1=wlntfp1;
[yinstrwith1,xinstrwith,meanny,meannx,meanty,meantx]=demean(yinstr1,xinstr,N,T,1);
results=ols(yinstrwith1,xinstrwith);
% Correction for degrees of freedom
results.sige=results.sige*(N*T-Nx)/(N*T-Nx-N-(T-1)); 
tmp = (results.sige)*(diag(results.xpxi));
results.tstat = results.beta./(sqrt(tmp));
% end correction
% vnames=strvcat('logp','pop','wpop','tax','wtax','comp','wcomp','logy','Wlogy');
vnames=strvcat('wlntfp','lnapp','lns','wlns','lnhr','wlnhr');
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


% Potential endogenous regressor 2
yinstr2=wlns;
[yinstrwith2,xinstrwith,meanny,meannx,meanty,meantx]=demean(yinstr2,xinstr,N,T,1);
results=ols(yinstrwith2,xinstrwith);
% Correction for degrees of freedom
results.sige=results.sige*(N*T-Nx)/(N*T-Nx-N-(T-1)); 
tmp = (results.sige)*(diag(results.xpxi));
results.tstat = results.beta./(sqrt(tmp));
% end correction
% vnames=strvcat('logp','pop','wpop','tax','wtax','comp','wcomp','logy','Wlogy');
vnames=strvcat('wlntfp','lnapp','lns','wlns','lnhr','wlnhr');
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
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,[x eta2resid],N,T,1);
K=3;

% 2SLS regression
fprintf(1,'2SLS model\n');
yinstrwith=[yinstrwith1 yinstrwith2];
results=tsls(ywith,yinstrwith,xwith,xinstrwith);
vnames=strvcat('lntfp','wlntfp','wlns','lnhr','lns');
prt_reg(results,vnames);
resiv=results.resid;
% Validity check IV
fprintf(1,'Validity IV\n');
results=ols(resiv,xinstrwith); % 
vnames=strvcat('resiv','lnapp','lns','wlns','lnhr','wlnhr');
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