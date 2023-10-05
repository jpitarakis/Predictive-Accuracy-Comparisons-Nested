
function [S0, S0_nw, S0_adj, S0_adj_nw, pv_S0, pv_S0_nw, pv_S0_adj, pv_S0_adj_nw] = Nested_Stats_S0(ehat1,ehat2,lam10,lam20)
    
% ==========================================================================================================================
% This function outputs the S(lam10,lam20) test statistics and their pvalues 
% as developed in "A Novel Approach to Predictive Accuracy Testing in Nested Environments" 
% Econometric Theory, 2023, Jean-Yves Pitarakis. 
%  
% Inputs: 
%
% ehat1: sequence of forecast errors from model 1 (small model)
% ehat2: sequence of forecast errors from model 2 (larger (nesting) model)
% lam10: user-input in (0,1] e.g., lam10 = 1
% lam20: user-input in (0,1], lam10 \neq lam20 e.g., lam20 = 0.90
%
% Outputs: 
% 
% S0: unadjusted test statistic and its pvalue pv_S0 (homoskedastic standardization)
% S0_nw : unadjusted test statistic and its pvalue pv_S0_nw (hac standardization)
% S0_adj: adjusted test statistic and its pvalue pv_S0_adj(homoskedastic standardization)
% S0_adj_nw : adjusted test statistic and its pvalue pv_S0_adj_nw (hac standardization)
% 
% 
% Remarks: (i) there are 4 variants of the same test statistic (unadjusted, 
% unadjusted+hac, adjusted, adjusted+hac). The hac versions of the statistics 
% use the Bartlett Kernel. The adjusted versions of the statistics use a modified 
% series of squared forecast errors "ehat2" that is size-neutral but power-enhancing 
% and is therefore preferred. (ii) all standard errors are based on the
% errors from the larger model (model 2). (iii) HAC standardizations use
% the COVNW function with Bartlett weights from Kevin Sheppards's MFE Matlab toolbox
% 
% ========================================================================

[n,~] = size(ehat1);
l10 = round(n*lam10);
l20 = round(n*lam20);

ehat1sq = ehat1.^2;
ehat2sq = ehat2.^2;
ehat2sq_adj = ehat2.^2-(ehat1-ehat2).^2;

mse_2 = mean(ehat2sq);
nuhat = ehat2sq-mse_2;

% Conditionally homoskedasticity based normalisers 

psisq = (nuhat)'*(nuhat)/n;
vcoef =(abs(lam10-lam20))/(lam10*lam20);
lrvar_chom = vcoef*psisq;
std_lrvar_chom = sqrt(lrvar_chom);


% Newey West based normalisers

%psisq_nw = lr_var(nuhat);
nlag = min(floor(1.2*n^(1/3)),n);
psisq_nw = covnw(nuhat,nlag,1);
lrvar_chet = vcoef*psisq_nw;
std_lrvar_chet = sqrt(lrvar_chet);


Z = (n/l10)*((sum(ehat1sq(1:l10))/sqrt(n))-(l10/l20)*(sum(ehat2sq(1:l20))/sqrt(n)));
Z_adj = (n/l10)*((sum(ehat1sq(1:l10))/sqrt(n))-(l10/l20)*(sum(ehat2sq_adj(1:l20))/sqrt(n)));

S0 = Z/std_lrvar_chom;
S0_nw = Z/std_lrvar_chet;

S0_adj = Z_adj/std_lrvar_chom;
S0_adj_nw = Z_adj/std_lrvar_chet;

pv_S0 = 1-normcdf(S0);
pv_S0_nw = 1-normcdf(S0_nw);

pv_S0_adj = 1-normcdf(S0_adj);
pv_S0_adj_nw = 1-normcdf(S0_adj_nw);

