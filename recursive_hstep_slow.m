function [zehat]=recursive_hstep_slow(y,X,pi0,h)

% Recursive forecasts and forecast errors using an expanding window approach 
%
% Inputs:
%
% y: nx1 outcome series
% X: nxp predictor matrix (the program adds an intercept)
% pi0: fraction of the sample to start recursions (e.g., pi0=0.5)
% h: forecast horizon
% 
% Recursively LS-fitted Model: y(t) = X(t-h) beta + u(t), t=h+1,...,n
%
% Outputs: 
% 
% sequence of h-steps ahead forecast errors (n-h-k0+1 x 1)

[n,p] = size(X);
k0 = round(n*pi0);
reg=[X,ones(n,1)];

bhat=nan(p+1,n-h-k0+1);
zehat = nan(n-h-k0+1,1);

for s=k0:(n-h)
bhat(:,s-k0+1)=reg(1:s-h,:)\y(h+1:s);
end

for i = k0:n-h
   zehat(i-k0+1) = y(i+h)-reg(i,:)*bhat(:,i-k0+1);
end
