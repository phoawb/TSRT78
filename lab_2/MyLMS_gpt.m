function [th,s_hat,err]=MyLMS_gpt(y,u,nb,nk,mu,lambda)
% [th,s_hat]=MyLMS(y,u,nb,nk,mu,lambda)
%
% A leaky LMS algorithm
%
% Inputs:       y       N x 1 vector with signal measurements (y=s+n)
%               u       N x 1 vector with input values
%               nb      model order
%               nk      number of samples the input should be delayed
%               mu      step length
%               lambda  leakage factor
%
% Outputs       th      N x (nb+1) matrix with the estimated filter parameters
%               s_hat   N x 1 vector with the filtered signal
%               err     N x 1 vector with the innovation, i.e., err(k)=y(k)-phi(k)'*theta(k-1)
%
% Author: YYY
% Date: XXX
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Place your code here
N = length(y);
th = zeros(N, nb+1);
s_hat = zeros(N, 1);
err = zeros(N, 1);

for k = 1:N
    if k > nk + nb
        phi = [1; zeros(nb, 1)]; % Initialize phi with zeros
        for j = 1:min(nb, k-nk)
            phi(j+1) = u(k-nk-j+1);
        end
        s_hat(k) = th(k-1, :) * phi; % Estimated output
        err(k) = y(k) - s_hat(k); % Error calculation
        th(k, :) = (1 - lambda) * th(k-1, :) + 2 * mu * err(k) * phi'; % Parameter update
    end
end
