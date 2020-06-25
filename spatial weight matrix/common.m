function [w]=common(x)

%% Generate spatial weight matrix of common group
%% the weight matrix is equivalent to within group average
%% x is the column of group index

[N,T]=size(x);
%% make sparse matrix of common group 
for i=1:N
for j=1:N
if j==i
d(i,j)=0;
elseif x(j)==x(i)
d(i,j)=1;
else
d(i,j)=0;
end
end
end
 
% normalize the row vectors
for i=1:N
    if sum(d(i,:),2)==0
        w(i,:)=d(i,:);
    else
        w(i,:)=d(i,:)/sum(d(i,:),2);
    end
end