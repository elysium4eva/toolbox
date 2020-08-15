function [w]=head(N,M)
%% Generate M ahead and M behind  matrix spatial weight matrix N unit
%% N is the number of spatial units
%% N>2M+1

d(1,:)=[0 ones(1,M) zeros(1,N-2*M-1) ones(1,M)];
for i=2:N
d(i,1)=d(i-1,N);
for j=2:N
d(i,j)=d(i-1,j-1);
end
end

% normalize the row vectors
for i=1:N
w(i,:)=d(i,:)/sum(d(i,:),2);
end