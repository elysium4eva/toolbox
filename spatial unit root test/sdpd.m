function results = sdpd(y,W,info)
%% estimate heterogeneous spatial dynamic model using QMLE
%% y_it=pi+rho*W*y_it+gamma_i*y_it-1+lambda_i*W*y_it-1+uit
%% info.model=0 without intercept and trend
%% info.model=1 with intercept only
%% info.model=2 with intercept and trend
%% result.rho is the concentratd MLE of rho
%% result.beta is the estimated beta
%% result.sigma is the esimated std error


[N,T]=size(y);
trend=2:T;


if info.model==0
%% ML estimate of theta
parmold=0;
options.Display='off';
options.MaxFunEvals=1000;
options.MaxIter=1000;
options.TolX=0.001;
options.TolFun=0.001;
[rho]=fminsearch('f_sdpd',parmold,options,y,W,info); 

k=3;
beta=zeros(N,2);
sigma=zeros(N,1);
Wy=W*y;

for i=1:N
    x=[y(i,1:T-1)' Wy(i,1:T-1)'];
    Mx=eye(T-1)-x*(x'*x)^(-1)*x';
    sy=y(i,2:T)-rho*Wy(i,2:T);
    beta(i,:)=(x'*x)^(-1)*x'*sy';
    sigma(i)=sy*Mx*sy'/(T-k);
end
results.rho=rho;
results.beta=beta;
results.sigma=sigma;

elseif info.model==1
        %% ML estimate of theta
parmold=0;
options.Display='off';
options.MaxFunEvals=1000;
options.MaxIter=1000;
options.TolX=0.001;
options.TolFun=0.001;
[rho]=fminsearch('f_sdpd',parmold,options,y,W,info); 

k=4;
beta=zeros(N,3);
sigma=zeros(N,1);
Wy=W*y;

for i=1:N
    x=[ones(T-1,1) y(i,1:T-1)' Wy(i,1:T-1)'];
    Mx=eye(T-1)-x*(x'*x)^(-1)*x';
    sy=y(i,2:T)-rho*Wy(i,2:T);
    beta(i,:)=(x'*x)^(-1)*x'*sy';
    sigma(i)=sy*Mx*sy'/(T-k);
end
results.rho=rho;
results.beta=beta;
results.sigma=sigma;

else
%% ML estimate of theta
parmold=0;
options.Display='off';
options.MaxFunEvals=1000;
options.MaxIter=1000;
options.TolX=0.001;
options.TolFun=0.001;
[rho]=fminsearch('f_sdpd',parmold,options,y,W,info); 

k=5;
beta=zeros(N,4);
sigma=zeros(N,1);
Wy=W*y;

for i=1:N
    x=[ones(T-1,1) trend' y(i,1:T-1)' Wy(i,1:T-1)'];
    Mx=eye(T-1)-x*(x'*x)^(-1)*x';
    sy=y(i,2:T)-rho*Wy(i,2:T);
    beta(i,:)=(x'*x)^(-1)*x'*sy';
    sigma(i)=sy*Mx*sy'/(T-k);
end
results.rho=rho;
results.beta=beta;
results.sigma=sigma;
end


    
        
end
