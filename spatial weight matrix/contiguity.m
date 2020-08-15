function [Wc] = contiguity(N,T,list)
%   Calculate the contiguity matrix from list of neighbor regions 
%   Not row normalized


list=list(1:2*N,:);
 for n=1:N
 m(n,:)=horzcat(list(2*n-1,:),list(2*n,:));
 end

m1=sortrows(m,1);
% Reform the matrix "list" into a matrix with the first column the index of region, second column the number of contiguous regrions
% And last columns are the indexes of contiguous regions of region i
ww(1,:)=zeros(1,N);

for i=1:N
    if m1(i,2)==0
       ww(i,:)=zeros(1,N);
    else
mi=m1(i,[T+1:2*T]);
y = 1:1:N;
for t=1:T
p=find(y==mi(t));
   ww(i,p)=1;
end
end
end


for i=1:N
    if sum(ww(i,:),2)==0
        for j=1:N
            if j==i
                w(i,j)=0;
            else
        w(i,j)=1/(N-1);
            end
        end
    else
        w(i,:)=ww(i,:);
    end
end
Wc=w;