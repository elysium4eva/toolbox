function RRho= mcsar1(K,rho,R,info)
a=[10,15,20,30,50,75,100,200,500];
table(1,:)=a;
table(2,:)=ones(1,9)*50;
table(3,:)=ones(1,9);

[nv,m]=size(table);

Rho=zeros(1,K);
RRho=zeros(R,m);
for r=1:R
for j=1:m
for k=1:K
    result=sdpdunitroot(table(2,j),2*r,table(1,j),rho,info);
    Rho(k)=result.rho;
end   
RRho(r,j)=mean(Rho);  
end
end





