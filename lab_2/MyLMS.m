function [th,s_hat,err]=MyLMS_2(y,u,nb,nk,mu,lambda)
% [th,s_hat]=MyLMS(y,u,nb,nk,mu,lambda)
%
% A leaky LMS algorithm
%
% Inputs:       y       N x 1 vector with signal measurements (y=s+n)
%               u       N x 1 vector with input values
%               nb      model order
%               nk      number of samples the input should be delayed
%               mu      step length
%               lambda  leakage factor (actually 'gamma' in the book)
%
% Outputs       th      N x (nb+1) matrix with the estimated filter parameters
%               s_hat   N x 1 vector with the filtered signal
%               err     N x 1 vector with the innovation, i.e., err(k)=y(k)-phi(k)'*theta(k-1)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Place your code here

