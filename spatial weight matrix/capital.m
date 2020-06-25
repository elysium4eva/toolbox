function [wcap] = capital(list)
%CAPITAL construct a matrix for a star network with capital city in its
%center for each country
%list is N by 2 matrix, N is the number of regions
%The first column of list is the city id, id=0, the capital, id<>0,
%otherwise, the second column of list is the country id

[N,k]=size(list);

for i=1:N
    if list(i,1)==0
        for j=1:N
     if list(j,2)==list(i,2)
         w(i,j)=1;
     else 
         w(i,j)=0;
     end
        end
    else
        for j=1:N
             if list(j,1)==0 && list(j,2)==list(i,2)
                 w(i,j)=1;
             else
                 w(i,j)=0;
             end
        end
    end
    
end

%% row normalize w
for i=1:N
wcap(i,:)=w(i,:)/sum(w(i,:),2);
end
end
