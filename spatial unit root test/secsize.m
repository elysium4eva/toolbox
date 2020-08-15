function [results] = secsize(K,N,theta,info)
% empirical size of secsf test
% K is the number of trials 
% N is the number of weight matrix sparseness
table=[30,30,30,50,50,50;50,100,200,50,100,200;-1.699,-1.735,-1.759,-1.645,-1.688,-1.706];
table1=[30,30,30,50,50,50;50,100,200,50,100,200;-1.736,-1.676,-1.607,-1.676,-1.616,-1.525];
table2=[30,30,30,50,50,50;50,100,200,50,100,200;-2.314,-1.956,-1.383,-2.264,-1.895,-1.328];
cips=zeros(K);
size=zeros(N,6);

if info.model==0

for t=1:6
for n=1:N
for k=1:K
results=secsf(table(1,t),table(2,t),n,theta);
cips(k)=results.cips;
end
count=sum(cips(:)<table(3,t));
size(n,t)=count/K;
end
end
results.size=size;


elseif info.model==1
        
for t=1:6
for n=1:N
for k=1:K
results=secsfintercept(table(1,t),table(2,t),n,theta);
cips(k)=results.cips;
end
count=sum(cips(:)<table1(3,t));
size(n,t)=count/K;
end
end
results.size=size;

else
for t=1:6
for n=1:N
for k=1:K
results=secsftrend(table(1,t),table(2,t),n,theta);
cips(k)=results.cips;
end
count=sum(cips(:)<table2(3,t));
size(n,t)=count/K;
end
end
results.size=size;       
end
        
end




