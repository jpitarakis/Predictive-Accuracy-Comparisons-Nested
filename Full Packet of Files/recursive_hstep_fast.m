function [ehat]=recursive_hstep_fast(y,X,pi0,h)

% Recursive LS forecasting using an expanding window approach. 
% 
% Inputs:
%
% y: nx1 outcome series
% X: nxp predictor matrix
% pi0: fraction of the sample to start recursions (e.g., pi0=0.5)
% h: forecast horizon
% 
% 
% Recursively LS-fitted Model (with intercept): y(t) = X(t-h) beta + u(t), t=h+1,...,n
%
% Notes: (1) first estimation window is [1,...,k0] and last window is 
% [1,....,n-h] for k0 = round(n*pi0). First forecast is yhat(k0+h|k0)
% and last forecast is yhat(n|n-h). There are a total of (n-h-k0+1)
% forecasts and corresponding forecast errors. (2) this fast version of the
% recursive least squares algorithm uses the Sherman-Morrison matrix
% formula to avoid matrix inversions at each recursion. 
%
% Outputs: 
%
% sequence of h-steps ahead forecast errors (n-h-k0+1 x 1)
%
%  

[n,p_x] = size(X);
p = p_x+1;
X = [X,ones(n,1)];
k0 = round(n*pi0);

% initial estimate using data t=1,...,k0

iXmat_k0 = inv(X(1:k0-h,:)'*X(1:k0-h,:));
bhat_k0 = iXmat_k0*(X(1:k0-h,:)'*y(h+1:k0));

M = zeros(p,p,n-k0-h+1); % initialisations for storage
bhat = zeros(p,n-h-k0+1); % initialisations for storage

bhat(:,1) = bhat_k0;
M(:,:,1) = iXmat_k0;

for t=k0:n-h
   M(:,:,t-k0+2) = M(:,:,t-k0+1)-(M(:,:,t-k0+1)*X(t+1-h,:)'*X(t+1-h,:)*M(:,:,t-k0+1)/(1+X(t+1-h,:)*M(:,:,t-k0+1)*X(t+1-h,:)'));
   bhat(:,t-k0+2) = bhat(:,t-k0+1)+M(:,:,t-k0+2)*X(t+1-h,:)'*(y(t+1)-X(t+1-h,:)*bhat(:,t-k0+1));
end

ehat = zeros(n-h-k0+1,1);

for s=k0:n-h
ehat(s-k0+1) = y(s+h)-X(s,:)*bhat(:,s-k0+1);
end



