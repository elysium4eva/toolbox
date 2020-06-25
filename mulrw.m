%% make n random walk series with innovation of normal distribution with standard error sigma for T periods
%% graph the series


function y=mulrw(mu,sigma,n,T)
x = linspace(1,T,T);   

for i=1:n
y(i,1)=0;
for t=1:T-1
    y(i,t+1)=normrnd(mu,sigma)+y(i,t);
end

plot(x,y(i,:),'-s');
hold on
end


