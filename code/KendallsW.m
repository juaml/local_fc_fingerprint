% taken from: https://github.com/ekmolloy/fmri_test-retest/blob/master/MATLAB/KendallsW.m
function [W, Chi2, df, pval] = KendallsW(X)
% ---------------------------------------------------------------
% Kendall's Correlation of Concordance (W)
%
% [W] = KendallsW(X) 
%       returns the Kendall's W of X
% [W, Chi2, df, pval] = KendallsW(X) 
%                       returns the Kendall's W of X with 
%                       Friedman's Chi-Squared, degrees of 
%                       freedom, and p-value
%
% The input X is an n x p matrix, where n is the number of 
% objects and p is the number of judges (i.e., ranking occurs 
% within columns).
% ---------------------------------------------------------------
% Script:
   [n,p] = size(X);

   poolsize = 0;
   pl = gcp('nocreate');
   if ~isempty(pl)
      poolsize = pl.NumWorkers;
   end
   
   % Tie Correction Factors
   T = nan(1,p);
   parfor(i = 1:p, poolsize)
        [V,iv,ix] = unique(X(:,i)); ixs = sort(ix);
        t = hist(ixs, ixs(end));
        T(i) = sum(t.^3 - t);
   end

   % Kendall's W
   R = sum(tiedrank(X), 2); 
   R_mean = mean(R);
   S = sum((R - R_mean).^2);
   W = (12*S) / ((p^2*(n^3-n)) - p*sum(T));

   % Friedman's Chi-Squared, degrees of freedom, and p-value
   Chi2 = p*(n-1)*W; df = n-1;
   pval = gammainc(Chi2/2, df/2, 'upper'); 
% ---------------------------------------------------------------
