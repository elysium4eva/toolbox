function plottop(d,K)
%PLOTTOP plot sum of top K panels for T
%   Detailed explanation goes here
[T,N]=size(d);
A(1)=0;
for t=1:T
    dt=sort(d(t,:),'descend');
    A(t)=sum(dt(1:K))/N;
end

plot(A);