function [results]=gravitynonnorm(w1,w2,w3,m_o,m_d,alpha,N,T)
%% The purpose of the function: using estimated alpha to contruct spatial
%% weight matrix of gravity form, the entry
%% w_ij=w1_ij^alpha(1)*w2_ij^alpha(2)*w3_ij^alpha(3)*m_0^alpha(4),...
%% w1 is matrix with dummy entries  
%% m_o is a N*T by k1 matrix of origin variables
%% m_d is a N*T by k2 matrix of destination variables
%% alpha is the parameter vector of number of constant weight matrices kw=3 plus number of variables in destination and origin

[Nob,k1]=size(m_o);
[Nob,k2]=size(m_d);
[k,junk]=size(alpha);

mo1=log(m_o(1:N,:));
md1=log(m_d(1:N,:));
for i=1:N
    for j=1:N
        if j==i 
            A(i,j)=0;
        else
            A(i,j)=exp(alpha(1)*log(w1(i,j))+alpha(2)*w2(i,j)+alpha(3)*w3(i,j)+alpha([4:3+k1])'*mo1(i,:)'+alpha([4+k1:3+k1+k2])'*md1(j,:)');
        end
    end
end 
lambda=eig(A);
A=A/max(lambda); 
d(1,:)=sum(A);

for t=2:T
mo1=log(m_o((t-1)*N+1:t*N,:));
md1=log(m_d((t-1)*N+1:t*N,:));
for i=1:N
    for j=1:N
        if j==i 
            B(i,j)=0;
        else
            B(i,j)=exp(alpha(1)*log(w1(i,j))+alpha(2)*w2(i,j)+alpha(3)*w3(i,j)+alpha([4:3+k1])'*mo1(i,:)'+alpha([4+k1:3+k1+k2])'*md1(j,:)');
        end
    end
end
lambda=eig(B);
B=B/max(lambda); 
A=blkdiag(A,B);
d(t,:)=sum(B);
end

results.w=A;
results.d=d;

