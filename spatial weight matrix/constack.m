function [w]=constack(N,T,list)
% Purpose: transfer the contiguity list of N regions into contiguity matrix
% "list" is the spreadsheet of contiguity list with odd rows the region's index and number of contiguous regions£¬even rows are the indexes of region's contiguous regions
% T is the maximum number of contiguous regions

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
w(i,:)=ww(i,:)/sum(ww(i,:),2);
    end
end

% Search in contiguous regions of region i, if there is a j, then assign w_ij=1