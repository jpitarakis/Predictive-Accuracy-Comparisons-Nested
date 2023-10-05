function [y,X] = dgp(n,p,k,sigsq_u,correl_uv,Omega_vv,Phimat,beta0,beta1)

% 
% This function simulated data from the following DGP: 
% ========================================================================
% y(t) = beta0 + x(t-1)'beta1 + u(t)
% x(t) = Phimat(p,p,1) x(t-1) + ... + Phimat(p,p,k) x(t-k) + v(t)
%
% Inputs:
%
% p = number of predictors in the predictive regression model (generated from a p-dimensional VAR(k) model). 
% n = sample size
% k = lag length of the VAR 
% sigsq_u = scalar variance of the u(t)
% Omega_vv = p x p vcov of the p-dimensional v(t)=(v_1(t),...,v_p(t))' vector i.e. E[v(t)v(t)']
% correl_uv = 1 x p vector of correlations between u(t) and (v_1(t),...,v_p(t)) i.e. E[u(t) v(t)] 
% Phimat: (p,p,k) matrices of coefficients of the VAR (just (p,p) if k=1)
% beta0: scalar intercept
% beta1: p x 1 vector of slope parameters (if set to a px1 vector of zeros the DGP becomes y(t) = b0 + u(t) 
%
% Outputs:
%
% Outcome series y. 
%
% Example usage: y(t) = 1 + 0.5 x1(t-1) + 0.25 x2(t-1) + u(t) 
% X(t) = (x1(t),x2(t)
% X(t) = Phimat X(t-1) + v(t), Phimat = [0.25,0.0;-0.50,0.10]
% Omega_vv = [1,0;0,1]
% omega_uv = [0.20 -0.20]
% sigsq_u = 1
% beta0 = 1
% beta1 = [0.5;0.25]
% n = 500
% 
% ================================================================

y = zeros(n,1);
omega_uv = correl_uv.*sqrt(sigsq_u).*(sqrt(diag(Omega_vv)))';
Sigma = [sigsq_u, omega_uv ; omega_uv', Omega_vv];
errors = mvnrnd(zeros(p+1,1), Sigma, n); 
u = errors(:,1);

X = varsim(p, Phimat, Omega_vv, k, n);
[nx,~] = size(X);
y(2:n) = beta0 + X(1:nx-1,:)*beta1+u(2:nx);
