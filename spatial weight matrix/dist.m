function [d]=dist(x,q)
%DIST Summary of this function goes here
%% This function calculate the inverse distance matrix without normalization
%% Generate spatial weight matrix of distance from coordinates of geographic unit
%% N is the number of units
%% The first column is the latitude, second column is the longitude

[N,T]=size(x);
R=6371;
X=x*pi/180;

%% make inverse distance matrix power by q
for i=1:N
for j=1:N
if j==i
d(i,j)=0;
else
d(i,j)=(acos(sin(X(i,1))*sin(X(j,1))+cos(X(i,1))*cos(X(j,1))*cos(X(j,2)-X(i,2)))*R)^(-q);
end
end
end
 
end

