function cv= mcsec(K,theta)
a=[10,15,20,30,50,70,100,200];
table(1,:)=reshape(ones(8,1)*a,[64,1])';
table(2,:)=reshape(a'*ones(1,8),[64,1])';
table(3,:)=0.2*table(2,:);

[nv,m]=size(table);
cv=zeros(3,m);

cips=zeros(1,K);
for j=1:m
for k=1:K
    result=secsfintercept(table(2,j),table(1,j),table(3,j),theta);
    cips(k)=result.cips;
end
cv(1,j)=prctile(cips,1);
cv(2,j)=prctile(cips,5);    
cv(3,j)=prctile(cips,10);    
end

