function TTheta= mcsec3(K,theta,R)
a=[5,10,15,20,30,50,75,100,200,500];
table(1,:)=a;
table(2,:)=ones(1,10)*50;
table(3,:)=ones(1,10);

[nv,m]=size(table);

Theta=zeros(1,K);
TTheta=zeros(R,m);
for r=1:R
for j=1:m
for k=1:K
    result=secsfintercept(table(2,j),table(1,j),2*r,theta);
    Theta(k)=result.thetahat;
end   
TTheta(r,j)=mean(Theta);  
end
end



