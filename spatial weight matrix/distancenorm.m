function [w]=distancenorm(x,q)
%% Generate spatial weight matrix of distance from coordinates of geographic unit
%% N is the number of units
%% The first column is the latitude, second column is the longitude
%% The distance is simpler euclidean distance of two regions, which allows repeated observation in the same location 
%% q is the power coefficient of inverse distance

[N,T]=size(x);
for i=1:N
for j=1:N
if j==i
d(i,j)=0;
else
d(i,j)=([x(i,1)-x(j,1)]^2+[x(i,2)-x(j,2)]^2+1)^(-0.5*q);
end
end
end

% normalize the row vectors
for i=1:N;
w(i,:)=d(i,:)/sum(d(i,:),2);
end