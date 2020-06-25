%%Illustration of FCLT/CMT

sigma=1; 
T=200;
N=1000;

for i=1:N
y(i,:)=randomwalk(sigma,T);
end

%% make a vector of trials of random walks at t=200
t=200;
y1=y(:,t)/T^0.5;

%% make a empirical probability density distribution for random trials of random walks
[f,xi] = ksdensity(y1);
plot(xi,f);
hold on



%% plot a normal distribution with standard deviation 
Mu=0;
Sigma=(t/T)^0.5;
x = -3:.1:3;
F= exp(-(x-Mu).^2./(2*Sigma^2))./(Sigma*sqrt(2*pi));
plot(x,F,'LineWidth',1);

hold on

