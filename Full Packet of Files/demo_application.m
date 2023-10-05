
% 
% This demo program demonstrates the use of the matlab codes accompanying 
% "A Novel Approach to Predictive Accuracy Testing in Nested Environments" 
% Econometric Theory, 2023, Jean-Yves Pitarakis. 
%
% The program loads a "dummy" data set from an excel file containig three
% time series: y, x1 and x2. 
%

[numbers, strings, raw] = xlsread('dummy_data_2.xlsx');
y = numbers(:,1);
x1 = numbers(:,2);
x2 = numbers(:,3);


% Predictive accuracy between Model 1 and Model 2
%
% Model 1: y(t) = b0 + b1 x1(t-1) + u(t)
% Model 2: y(t) = b0 + b1 x1(t-1) + b2 x2(t-1) + u(t)
%
% Obtaining the recursive forecast errors starting from a initial fraction
% of 25% of the observations (i.e. pi0 = 0.25)
% 
%

pi0 = 0.25;
ehat1 = recursive_hstep_fast(y,[x1],pi0,1);
ehat2 = recursive_hstep_fast(y,[x1 x2],pi0,1);

lam10o = 1; lam20o = 0.9;
[S0, S0_nw, S0_adj, S0_adj_nw, pv_S0, pv_S0_nw, pv_S0_adj, pv_S0_adj_nw] = Nested_Stats_S0(ehat1,ehat2,lam10o,lam20o);
% S0(lam1,lam2) test statistics and their pvalues 
out_S0 = [S0, S0_nw, S0_adj, S0_adj_nw; pv_S0, pv_S0_nw, pv_S0_adj, pv_S0_adj_nw]


lam20 = 1; tau0 = 0.8;
[Sbar, Sbar_nw, Sbar_adj, Sbar_adj_nw, pv_Sbar, pv_Sbar_nw, pv_Sbar_adj, pv_Sbar_adj_nw] = Nested_Stats_Sbar(ehat1,ehat2,lam20,tau0);
% Sbar(tau0;lam20) test statistics and their pvalues 
out_Sbar = [Sbar, Sbar_nw, Sbar_adj, Sbar_adj_nw; pv_Sbar, pv_Sbar_nw, pv_Sbar_adj, pv_Sbar_adj_nw]



