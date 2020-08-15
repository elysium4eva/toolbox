function results= mcsec2(K,theta)
a=[6,7,8,9,10,15,20,25,30,40,50,100,200,500,1000];
table(1,:)=a;
table(2,:)=ones(1,15)*50;
table(3,:)=0.2*table(2,:);

[nv,m]=size(table);
results.sta=zeros(3,m);

cips=zeros(1,K);
Theta=zeros(1,K);
for j=1:m
for k=1:K
    result=secsfintercept(table(2,j),table(1,j),table(3,j),theta);
    cips(k)=result.cips;
    Theta(k)=result.thetahat;
end
results.sta(1,j)=mean(cips);
results.sta(2,j)=std(cips);    
results.sta(3,j)=mean(Theta);    
end


