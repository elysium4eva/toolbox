function [results] = gravitycons(w1,w2,w3,alpha,N,T)
%Calculate gravity weight matrix using only time invariant matrix
%Detailed explanation goes here
for i=1:N
    for j=1:N
        if j==i 
            A(i,j)=0;
        else
            A(i,j)=exp(alpha(1)*log(w1(i,j))+alpha(2)*w2(i,j)+alpha(3)*w3(i,j));
        end
    end
A(i,:)=A(i,:)/sum(A(i,:));
end
d(1,:)=sum(A);

for t=2:T
for i=1:N
    for j=1:N
        if j==i 
            B(i,j)=0;
        else
            B(i,j)=exp(alpha(1)*log(w1(i,j))+alpha(2)*w2(i,j)+alpha(3)*w3(i,j));
        end
    end
    B(i,:)=B(i,:)/sum(B(i,:));
end
A=blkdiag(A,B);
d(t,:)=sum(B);
end

results.w=A;
results.d=d;
end

