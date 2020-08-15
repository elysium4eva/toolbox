function results = est_delta(dTN,M)
% PURPOSE: estimates the degrees of dominance for the top M units
% INPUTS: dTN = a T by N matrix of outdegrees
%         M = number of top pervasive units to be reported 
% OUTPUT: results = a structure, where
%         results.deltaM = degrees of dominance for the top M units
%         results.idx_deltaM = index array corresponding to the estimates
%         results.Nnz = number of units with non-zero outdegrees
%         results.se_deltaM = standard error of the estimates 
%                (Note that the standard error cannot be computed when T=1.)
% ------------------------
% NOTE: This code was written for large N and short T panels (including T=1)
%       under the assumption that the errors in the exponent specification 
%       are IID. The extension to large N and T panels was considered in an
%       Online Supplement.
%-------------------------
% Please let me know if you have any questions using this code. All
% feedback is much appreciated.
% Cynthia Fan Yang
% Email: yang488@usc.edu
% Last updated: August 2017

T = size(dTN,1);
lndTN = log(dTN);
lndTN(isinf(lndTN)) = 0;  % replace Inf with 0 
Tnz = sum(lndTN~=0,1);  % count nonzero Ts for each i 
Tnz = Tnz + (Tnz==0);   % protect against zeros
lndi_bar = sum(lndTN,1)./Tnz;
Nnz = sum(lndi_bar~=0);
lnd_bar = sum(lndi_bar)/Nnz;
deltai_hat = (lndi_bar-lnd_bar)./log(Nnz); 
[deltai_sorted, idx_deltai] = sort(deltai_hat,'descend'); 
deltaM = deltai_sorted(1:M);
idx_deltaM = idx_deltai(1:M);
results.deltaM = deltaM;        
results.idx_deltaM = idx_deltaM; 
results.Nnz = Nnz;          

% Compute standard error if T>1
if T==1
    results.se_deltaM = 'N/A';
elseif T>1
    vhatTN = bsxfun(@minus, lndTN - lnd_bar, deltai_hat.*log(Nnz));
    if min(Tnz)==T   % balanced panel
        sig2v = sum(sum(vhatTN.^2))/(T-1)/Nnz;   
        var_deltaM = sig2v/log(Nnz)^2/T*(1-1/Nnz);
    else             % unbalanced panel
        vhatTN = bsxfun(@minus, lndTN - lnd_bar, deltai_hat.*log(Nnz));
        vhatTN = vhatTN.*(lndTN~=0);  
        idx_Tnz1 = (Tnz==1);
        vhatTN_copy = vhatTN;
        Tnz_copy = Tnz;
        Tnz_copy(idx_Tnz1) = [];
        vhatTN_copy(:,idx_Tnz1) = [];
        sig2v = sum(sum(vhatTN_copy.^2)./(Tnz_copy-1))/size(Tnz_copy,2); 
        Tdomi = Tnz(idx_deltaM(1));    % sample size of the most dominant unit
        var_deltaM = sig2v/log(Nnz)^2/Tdomi*(1-1/Nnz);
    end
    se_deltaM = sqrt(var_deltaM);
    results.se_deltaM = se_deltaM;
end

end