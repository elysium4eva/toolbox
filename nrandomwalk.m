function y=nrandomwalk(sigma,T,N) 
%% plot N random walks to show the non-stationarity
for i=1:N
    y(i,:)=randomwalk(sigma,T);
    plot(y(i,:)');
    hold on
end