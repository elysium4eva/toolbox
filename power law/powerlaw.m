function results=powerlaw(d,M)

%% The purpose of this function is to report estimation of power law coefficients given a matrix of network degree d by T periods of observations
%% d is an N by T matrix of network outdegree, pct_set is the cut-off value for the first two estimators
%% beta_GI reports Gabaix and Ibragimov (2011) estimator by OLS, with its standard error sigmaGI  
%% beta_MLE reports Hill estimator (Hill et al., 1975), as the ML estimator, standard error sigmaMLE
%% beta_CSN reports the feasible ML estimator of Clauset et al. (2009, CSN)
%% beta_Max reports the inverse of the extremum estimator delta (Yang 2020)
%% M = number of top pervasive units to be reported 

[N,T]=size(d);
pct_set = [0.1,0.2,0.3];   % assuming cut-off values: 10%, 20%, 30%
len_pct = length(pct_set);


%========== Extremum estimator of delta ===========
results = est_delta(d',M);
    
for t=1:T

for pi= 1:3  % under the assumed cut-off value
        pct_cut = pct_set(pi);      
        dN_sorted = sort(d(:,t),'descend');
        cut = round(pct_cut*N);
% GI regresion -- Log-log regression with Gabaix-Ibragimov correction
        b_GI = regress(log((1-0.5:1:cut-0.5)'), [ones(cut,1),log(dN_sorted(1:cut))]);
        [se_b_GI] = (sqrt(2/cut))*b_GI(2);
        results.beta_GI(t,pi) = -b_GI(2);
        results.se_beta_GI(t,pi) = -se_b_GI;
% MLE
        results.beta_mle(t,pi) = cut/(sum(log(dN_sorted(1:cut))) - cut*log(dN_sorted(cut)));   
        results.se_beta_mle(t,pi) =results.beta_mle(t)/sqrt(cut);          
    end   

% Infeasible MLE (CSN estimation), calls function plfit.m by Clauset, Shalizi and Newman (2009)

    [alpha_CSN(t), xmin_CSN(t)] = plfit(nonzeros(dN_sorted),'nosmall');   
    ntail_CSN(t) = length(find(dN_sorted > xmin_CSN(t)));
    results.se_alpha_CSN(t) = (alpha_CSN(t)-1)/sqrt(ntail_CSN(t));
    results.beta_CSN(t)= alpha_CSN(t)- 1;
    results.pct_ntail_CSN(t) = ntail_CSN(t)/N;   % estimated cut-off value  
  
end


