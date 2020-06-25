function [wg]=gravity(coordinates,mass,q)
%% Generate spatial weight matrix of distance from coordinates of geographic unit
%% N is the number of units
%% The first column of coordinates is the latitude, second column is the longitude
%% q>0 is the power coefficient of inverse distance
%% mass is a panel variable of N spatial units and T time

[n,u]=size(mass);
[N,v]=size(coordinates);
T=n/N;
y=reshape(mass,[N,T]);
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

for i=1:N
dd(i,:)=d(i,:)/sum(d(i,:),2);
end


w=diag(y(:,1))*dd*diag(y(:,1));
for i=1:N
wg(i,:)=w(i,:)/sum(w(i,:),2);
end


% block diagnalize the gravity matrix for T period, and normalize the row vectors
for t=2:T
ww=diag(y(:,t))*dd*diag(y(:,t));
for i=1:N
www(i,:)=ww(i,:)/sum(ww(i,:),2);
end
wg=blkdiag(wg,www);
t=t+1;        
end

