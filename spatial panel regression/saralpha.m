function [results] = saralpha(y,x,w,N,T,info)
%Estimate the SAR panel model with distance decay parameter alpha
% w is the inverse distance matrix
 
alphaold=1;
iter=0;converge=1.0;criteria=0.0001;itermax=25;
options.Display='off';
options.MaxFunEvals=1000;
options.MaxIter=1000;
options.TolX=0.000001;
options.TolFun=0.000001;

[alpharesult,fval,exitflag]=fminsearch('f_sarfe_alpha',alphaold,options,y,x,w,T,info);

walpha=w.^alpharesult;
for i=1:N
Walpha(i,:)=walpha(i,:)/sum(walpha(i,:),2);
end
results.alpha=alpharesult;
results=sar_panel_FE(y,x,Walpha,T,info);
end