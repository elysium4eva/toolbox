function [results] =slxendo(y,x,eta,N,T,W)
%ENDOGENOUSWEIGHT Summary of this function goes here
% 2SLS  Estimate the spatial model with endogenous weight matrix
% y is the Nob by 1 vector dependent variable stacked by region 1, region 2,...,region N at
% period 1; region 1, region 2,..., region N at period 2;...
% w is the prespecified N by N weight matrix containing endougeous variable
% eta
% info.model = 0 pooled model without fixed effects (default, x may contain an intercept)
%            = 1 spatial fixed effects (x may not contain an intercept)
%            = 2 time period fixed effects (x may not contain an intercept)
%            = 3 spatial and time period fixed effects (x may not contain an intercept)

[nobs,nvar]=size(x); 

%% define the feasible instruments H
H=[W*x(:,1),eta,W*eta];
Z=[x];

%% define the projection matrix PH
PH=H*(H'*H)^(-1)*H';

JTbar=ones(T)/T;
ET=eye(T)-JTbar;
P=kron(JTbar,eye(N));
Q=eye(N*T)-P;
Hfe=Q*H;
%% define the projection matrix PH
PHfe=Hfe*(Hfe'*Hfe)^(-1)*Hfe';
%% delta is the parameter vector of 2sls estimator
results.beta=(Z'*PHfe*Z)^(-1)*Z'*PHfe*y;
results.resid=y-Z*results.beta;
results.sigu = results.resid'*results.resid;
results.sige = results.sigu/(nobs-nvar);



end


