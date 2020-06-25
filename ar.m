%% make a first order autoregressive series with innovation of normal distribution with standard error sigma for T periods
%% graph the series
%% the series is stationary if -1<rho<1

function y=rwdrift(mu,sigma,rho,T)
y(1)=0;
for t=1:T-1
    y(t+1)=normrnd(mu,sigma)+rho*y(t);
end
 x = linspace(1,T,T);   
 plot(x,y,'-s');