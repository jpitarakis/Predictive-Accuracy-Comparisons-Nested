
function [Sbar, Sbar_nw, Sbar_adj, Sbar_adj_nw, pv_Sbar, pv_Sbar_nw, pv_Sbar_adj, pv_Sbar_adj_nw] = Nested_Stats_Sbar(ehat1,ehat2,lam20,tau0)

% ==========================================================================================================================
% This function outputs the Sbar(tau0;lam20) test statistics and their pvalues 
% as developed in "A Novel Approach to Predictive Accuracy Testing in Nested Environments" 
% Econometric Theory, 2023, Jean-Yves Pitarakis. 
%  
% Inputs: 
%
% ehat1: sequence of forecast errors from model 1 (small model)
% ehat2: sequence of forecast errors from model 2 (larger (nesting) model)
% lam20: user-input in [0,1] e.g., lam20 = 1
% tau0 : user-input in (0,1] e.g., tau0 = 0.8
%
% Outputs: 
% 
% Sbar: unadjusted test statistic and its pvalue pv_Sbar (homoskedastic standardization)
% Sbar_nw : unadjusted test statistic and its pvalue pv_Sbar_nw (hac standardization)
% Sbar_adj: adjusted test statistic and its pvalue pv_Sbar_adj(homoskedastic standardization)
% Sbar_adj_nw : adjusted test statistic and its pvalue pv_Sbar_adj_nw (hac standardization)
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
l20 = round(n*lam20);
taua = round(n*tau0);
taub = round(n*(1-tau0));

ehat1sq = ehat1.^2;
ehat2sq = ehat2.^2;
ehat2sq_adj = ehat2.^2-(ehat1-ehat2).^2;

mse_2 = mean(ehat2sq);
nuhat = ehat2sq-mse_2;

psisq = nuhat'*nuhat/n;

if tau0 == 0
    vartau=(1+2*(lam20)*log(lam20))/(lam20);
elseif tau0>0 && tau0<1 && lam20<=tau0
    vartau = (((1-tau0)^2)+2*lam20*(1-tau0+log(tau0)))/(lam20*((1-tau0)^2));
 elseif tau0>0 && tau0<1 && lam20>tau0
     vartau = (1-(tau0^2)+2*lam20*((1-tau0)*log(lam20)+tau0*log(tau0)))/(lam20*((1-tau0)^2));
end

% Conditional homoskedastic normaliser

lrvar_chom = vartau*psisq;
std_lrvar_chom = sqrt(lrvar_chom);

% Newey-West normaliser

%psisq_nw = lr_var(nuhat);
nlag = min(floor(1.2*n^(1/3)),n);
psisq_nw = covnw(nuhat,nlag,1);
lrvar_chet = vartau*psisq_nw;
std_lrvar_chet = sqrt(lrvar_chet);

zvec = zeros(n,1);
zvec_adj = zeros(n,1);

for l1 = (taua+1):1:n
zvec(l1) = (n/l1)*((sum(ehat1sq(1:l1))/sqrt(n))-(l1/l20)*(sum(ehat2sq(1:l20))/sqrt(n)));
zvec_adj(l1) = (n/l1)*((sum(ehat1sq(1:l1))/sqrt(n))-(l1/l20)*(sum(ehat2sq_adj(1:l20))/sqrt(n)));
end

Zbar = sum(zvec)/taub;
Zbar_adj = sum(zvec_adj)/taub;

Sbar = Zbar/std_lrvar_chom;
Sbar_adj = Zbar_adj/std_lrvar_chom;
pv_Sbar = 1-normcdf(Sbar);
pv_Sbar_adj = 1-normcdf(Sbar_adj);

Sbar_nw = Zbar/std_lrvar_chet;
Sbar_adj_nw = Zbar_adj/std_lrvar_chet;
pv_Sbar_nw = 1-normcdf(Sbar_nw);
pv_Sbar_adj_nw = 1-normcdf(Sbar_nw);


