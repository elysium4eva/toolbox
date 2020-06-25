function results=granular(W,N1,N2);
%% the purpose of this function is to test the granularity condition of a prespecified spatial weight matrx based on the definition of Pesaran(2009),
%% Pesaran and Tosetti, 2007)
%% the pre-specified spatial weight matrix has degree N, when norm of W is in the asymptotic order of N^-1/2
 
for i=N1:N2
    wi=head();
results.a(i)=trace(wi'*wi)^0.5/i^(-0.5);
end
plot results.a
