function [w]=distancenband(x,d)
%% Generate spatial weight matrix of distance from coordinates of geographic unit
%% N is the number of units
%% The first column is the latitude, second column is the longitude
%% d>0 is the radius by km of defined neibourhood if d_ij <=d, then w_ij=1, otherwise w_ij=0



[N,T]=size(x);
R=6371;
X=x*pi/180;

%% make distance matrix l
for i=1:N
for j=1:N
if j==i
l(i,j)=0;
else
l(i,j)=acos(sin(X(i,1))*sin(X(j,1))+cos(X(i,1))*cos(X(j,1))*cos(X(j,2)-X(i,2)))*R;
end
end
end

%make sparse matrix c within radius d
for i=1:N
for j=1:N
if j==i
c(i,j)=0;
elseif l(i,j)<=d
c(i,j)=1;
else
c(i,j)=0;
end
end
end
 
% normalize the row vectors
for i=1:N
    if sum(c(i,:),2)==0
        w(i,:)=c(i,:);
    else
        w(i,:)=c(i,:)/sum(c(i,:),2);
    end
end