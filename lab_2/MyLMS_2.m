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
    % Length of the input signal
    N = length(y);

    % Matrix with theta, we have nb+1 coefficients
    % for every time step
    th = zeros(nb+1, N+1);

    % The estimated signal
    s_hat = zeros(N, 1);

    % The difference between the signal and the
    % estimated signal
    err = zeros(N, 1);
    
    % Loop through every sample in the signal
    for t = 1:N
        % Construct the input vector phi if the delay is not
        % greater than the current time index
        if(t-nk > 0)
            % Create phi, the regression vector 
            % from the input vector 'u' and
            % take 'nb' previous samples from 'u'
            % with a delay of 'nk'
            phi = [
                u(t-nk:-1:max(t-nk-nb, 1)); 
                zeros(nk+nb+1-t, 1)
                ]
        else
            phi = zeros(nb+1, 1);
        end
        % Update the filter coefficients of theta for the
        % current time step based on leakage from previous
        % time step and update in gradient direction
        th(:, t+1) = (1-lambda)*th(:, t) + mu*phi*(y(t) - phi'*th(:, t));
        % the error is y minus estimated signal
        err(t) = y(t) - phi'*th(:, t);
        % subtract the estimated noise
        s_hat(t) = y(t) - phi'*th(:, t+1); % shat = y - y_hat
    end
    % Removes the first column from the th matrix 
    % (which was used for initialization) and 
    % transposes the matrix so that each row 
    % corresponds to the filter coefficients at 
    % a specific time step
    th = th(:, 2:end)';
