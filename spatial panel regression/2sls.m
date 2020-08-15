function [results] = 2sls(y,x,w,info)
%2SLS Summary of this function goes here
% The purpose of this function is to estimate panel spatial autoregressive
% model by 2 stage least square mothed
% y is the dependent variable stacked by region 1, region 2,...,region N at
% period 1; region 1, region 2,..., region N at period 2;...
% w is the prespecified N by N weight matrix 
% The notations are taken from Baltagi 2013, p325
% info.model = 0 pooled model without fixed effects (default, x may contain an intercept)
%            = 1 spatial fixed effects (x may not contain an intercept)
%            = 2 time period fixed effects (x may not contain an intercept)
%            = 3 spatial and time period fixed effects (x may not contain an intercept)


[N,N]=size(w);
T=size(y,1)/N;
W=kron(eye(T),w);
Zmiu=kron(ones(T),eye(N));

%% define the feasible instruments H
H=[x,W*x,W*W*x];
Z=[x,W*y];

%% define the projection matrix PH
PH=H*(H'*H)^(-1)*H';

if info.model=0
%% delta is the parameter vector of 2sls estimator
result.beta=(Z'*PH*Z)^(-1)*Z'*PH*y;

else if info.model=1
JTbar=ones(T)/T;
ET=eye(T)-JTbar;
P=kron(JTbar,eye(N));
Q=eye(N*T)-P;
Hfe=Q*H;
%% define the projection matrix PH
PHfe=Hfe*(Hfe'*Hfe)^(-1)*Hfe';
%% delta is the parameter vector of 2sls estimator
results.beta=(Z'*PHfe*Z)^(-1)*Z'*PHfe*y;
results.resid=y-Z*result.beta;
results.sigu = results.resid'*results.resid;
results.sige = results.sigu/(nobs-nvar);
sige=results.sige*((N*T-K)/(N*T));
loglik=-N*T/2*log(2*pi*sige)-1/(2*sige)*results.resid'*results.resid
else if info.model=2
        Hbe=P*H;
        %% define the projection matrix PH
        PHbe=Hbe*(Hbe'*Hbe)^(-1)*Hbe';
%% delta is the parameter vector of 2sls estimator
result.beta=(Z'*PHbe*Z)^(-1)*Z'*PHbe*y;

end

