function [wgv] = gravityv(coordinates,y,q)
 
%% Generate spatial weight matrix of distance from coordinates of geographic unit
%% N is the number of units
%% The first column of coordinates is the latitude, second column is the longitude
%% q>0 is the power coefficient of inverse distance
%% mass is a vector N spatial units 

[N,u]=size(y); 

R=6371;
X=coordinates*pi/180;

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

ww=diag(y)*d*diag(y);
for i=1:N
wgv(i,:)=ww(i,:)/sum(ww(i,:),2);
end

end