function cv= mcsec1(K,theta)
a=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];

table(1,:)=ones(1,20)*50;
table(2,:)=ones(1,20)*100;
table(3,:)=a;

[nv,m]=size(table);
cv=zeros(3,m);

cips=zeros(1,K);
for j=1:m
for k=1:K
    result=secsf(table(2,j),table(1,j),table(3,j),theta);
    cips(k)=result.cips;
end
cv(1,j)=prctile(cips,1);
cv(2,j)=prctile(cips,5);    
cv(3,j)=prctile(cips,10);    
end


