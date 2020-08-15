
%% define the feasible instruments H
H=[x,W*x,W*W*x];
Z=[x,W*y];

%% define the projection matrix PH
PH=H*(H'*H)^(-1)*H';

JTbar=ones(T)/T;
ET=eye(T)-JTbar;
P=kron(JTbar,eye(N));
Q=eye(N*T)-P;
Hfe=Q*H;
%% define the projection matrix PH
PHfe=Hfe*(Hfe'*Hfe)^(-1)*Hfe';


W=W1;
%% delta is the parameter vector of 2sls estimator
results.beta=(Z'*PHfe*Z)^(-1)*Z'*PHfe*y;
results.resid=y-Z*results.beta;
results.sigu = results.resid'*results.resid;
results.sige = results.sigu/(nobs-nvar);
sige=results.sige*((N*T-nvar)/(N*T));
results.loglik=-N*T/2*log(2*pi*results.sige)-1/(2*results.sige)*results.resid'*results.resid;
