function W=rooksq(n,m);
% construct a contiguity list of a rectangle lattice of n by m;
% w is a spreadsheet of n times m rows list with the first column the index
% of the region, the index of contiguous regions in the following columns
% the total column is 5
N=n*m;
 for i=2:N
     w(i,1)=i;
     w(i,2)=i+1;
     w(i,3)=i-1;
     w(i,4)=i-m;
     w(i,5)=i+m;
 end
 for i=2:m-1
      w(i,1)=i;
      w(i,2)=i+1;
      w(i,3)=i-1;
      w(i,4)=0;
      w(i,5)=i+m;
 end
 for i=(n-1)*m+2:N-1
      w(i,1)=i;
      w(i,2)=i+1;
      w(i,3)=i-1;
      w(i,4)=i-m;
      w(i,5)=0;
 end
 for i=2:n-1
     w(i*(m-1)+1,1)=i*(m-1)+1;
     w(i*(m-1)+1,2)=i*(m-1)+2;
     w(i*(m-1)+1,3)=0;
     w(i*(m-1)+1,4)=i*(m-1)+1-m;
     w(i*(m-1)+1,5)=i*(m-1)+1+m;
 end
 for i=2:n-1
     w(i*m,1)=i*m;
     w(i*m,2)=0;
     w(i*m,3)=i*m-1;
     w(i*m,4)=(i-1)*m;
     w(i*m,5)=(i+1)*m;
 end
      w(1,1)=1;
      w(1,2)=2;
      w(1,3)=m+1;
      w(1,4)=0;
      w(1,5)=0;
      
      w(m,1)=m;
      w(m,2)=m-1;
      w(m,3)=2*m;
      w(m,4)=0;
      w(m,5)=0;
      
      w((n-1)*m+1,1)=(n-1)*m+1;
      w((n-1)*m+1,2)=(n-2)*m+1;
      w((n-1)*m+1,3)=(n-1)*m+2;
      w((n-1)*m+1,4)=0;
      w((n-1)*m+1,5)=0;
      
      w(N,1)=N;
      w(N,2)=(n-1)*m;
      w(N,3)=N-1;
      w(N,4)=0;
      w(N,5)=0;

% expend the contiguity list into contiguous weight matrix

v=w(:,1);
for i=1:N
for j=2:5
p=find(v==w(i,j)); 
W(i,p)=1;
end 
end
% normalize the row vectors

for i=1:N
W(i,:)=W(i,:)/sum(W(i,:),2);
end
