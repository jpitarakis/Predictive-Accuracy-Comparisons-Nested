function data = varsim(p, A, Sigma, k, n)

% This function simulates data from a p-dimensional VAR with k lags
%
% X(t) = A1 X(t-1)+...+ Ak X(t-k) + v(t) for v(t)=(v1(t),...,vp(t))'.
%
% Inputs: 
%
% p: the dimension of the process
% A: the (p by p by k) array of coefficients
% Sigma: the p by p covariance matrix of the errors E[v(v)v(t)']
% k: the lag order of the process
% n: the number of observations to generate
%
% Output: 
%
% data: the n by p matrix of generated data
%
% 
% Remarks: (i) For a VAR(k=1), the function expects A to be a (p x p)
% coefficient matrix. For higher lag orders A must be specified as a three
% dimensional array (p x p x k). (ii) Note that the output matrix 
% is (n x p) dimensional and not (n-k x p) as it subsumes the initial
% values. One may wish to replace data with data(k+1:end,.) by
% de-commenting the last line in this code. 
%
% Example: Simulate a 2-dimensional VAR(1) with coefficient matrix 
% A = [0.5,0.1;0.25,0.0] and identity covariance of the errors. 
% Assume 100 observations
% varsim(2, A, eye(2), 1, 100)
%
% ==================================================================

% Initialize the data matrix with zeros
data = zeros(n, p);

% Generate the errors from a multivariate normal distribution
errors = mvnrnd(zeros(1, p), Sigma, n);

% Generate the first k observations from the errors
data(1:k, :) = errors(1:k, :);

% Loop over the remaining observations
for i = k+1:n
    % Initialize a vector to store the linear combination of previous observations
    x = zeros(p, 1);
    
    % Loop over the lags
    for j = 1:k
        % Add the product of the coefficient matrix and the previous observation at lag j
        x = x + A(:, :, j) * data(i-j, :)';
    end
    
    % Add the error term to the linear combination
    data(i, :) = x' + errors(i, :);
    
end
 %data = data(k+1:end, :);

end 