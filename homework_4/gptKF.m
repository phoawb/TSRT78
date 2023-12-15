clear; close all; clc;
load Aug11.mat % Load the data file

meas = ykf; % get measurement
n_meas = length(meas); % sequence length
time = (1:n_meas)'; % time vector (column vector)

% State-space model setup
% x[t] = [a_t; b_t]
% Since a_t and b_t are time-varying, we consider them as random walk processes
% x[t+1] = F*x[t] + w[t], where w[t] is the process noise with covariance Q
F = [1, 0; 0, 1]; % State transition matrix (identity matrix for random walk)
Q = diag([0.01, 0.01]); % Process noise covariance (tune these values as needed)
R = 1 %var(meas - mean(meas)); % Measurement noise covariance (can be tuned)
phi = [1, time(1)]; % Measurement matrix (initially for t=1)

% Initial state and covariance
x_est = [0; 0]; % Initial state estimate (can be tuned)
P = 100 * eye(2); % Initial error covariance (large value to represent high initial uncertainty)

% Kalman filter implementation
x_estimated = zeros(2, n_meas); % To store estimated states
for t = 1:n_meas
    % Prediction
    P_pred = F * P * F' + Q;
    
    % Update the measurement matrix for current time
    phi = [1, time(t)]; % Measurement matrix for time t
    
    % Update
    y_pred = phi * x_est; % Predicted measurement
    e = meas(t) - y_pred; % Measurement residual (innovation)
    S = phi * P_pred * phi' + R; % Residual covariance
    K = P_pred * phi' * (S \ eye(size(S))); % Kalman gain
    x_est = x_est + K * e; % Updated state estimate
    P = (eye(2) - K * phi) * P_pred; % Updated error covariance
    
    % Store the estimated state
    x_estimated(:, t) = x_est;
end

% Extract estimated parameters
a_est = x_estimated(1, :);
b_est = x_estimated(2, :);
s_est = a_est + b_est .* time'; % Estimated signal

% Plotting

figure;
plot(time, x_estimated);
title("parameters")
legend("a_t", "b_t")

figure;
plot(time, meas, 'linewidth', 2); hold on; % plot measurement
plot(time, s_est, 'r', 'linewidth', 2); % plot estimated signal
legend('Measurements', 'Estimated Signal');
title('Signal Estimation using Kalman Filter');
xlabel('Time');
ylabel('Signal Value');


figure;
subplot(311); 
plot(time, a_est);
subplot(312); 
plot(time, b_est);
subplot(313);
plot(time, b_est .* time');

a2 = 0.9101 * ones(1, 300);

