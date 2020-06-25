%% make a random walk series with innovation of normal distribution with standard error sigma for T periods
%% graph the series

function y=rwdrift(mu,sigma,T)
y(1)=0;
for t=1:T-1
    y(t+1)=normrnd(mu,sigma)+y(t);
end
 x = linspace(1,T,T);   
 plot(x,y,'-s');