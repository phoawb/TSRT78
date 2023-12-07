rng(0)                 % For reproducibility

% Given parameters
A_true = 1 / sqrt(2);  % True amplitude
xi_true = 0.1;         % Angular frequency
phi_true = pi / 2;     % Phase shift
lambda_var = 0.25;     % Variance of the noise
N = 1000;              % Number of samples

% Generating the true signal
t = 0:(N-1);           % Discrete time steps
s = A_true * sin(xi_true * t + phi_true);

% Generating the observed signal with noise               
e =  sqrt(lambda_var).*randn(1,N);
y = s + e;

% Plot zone
figure(1)
hold on
plot(t,y)
plot(t,s,'LineWidth',2)
legend('Observed signal', 'True signal')
title('True signal and observation')
xlabel('sample step t')
print('./homework_2_comparison_1', '-dpng')

% ----- Least Squares Estimation ------

% Datamatris
A = [sin(t*xi_true); 
     cos(t*xi_true)];

% Get the x_1 and x_2
x = pinv(A')*y';

% Amplitude from the x
A_hat = sqrt(x(1).^2 + x(2).^2)

% Phase-shift from the x
phi_hat = atan2(x(2), x(1))
 
% Estimated signal
y_hat = A_hat * sin(xi_true * t + phi_hat);

total_error = sum((y_hat - s).^2)

figure(2)
subplot(3,1,1)
plot(t,s)
title('True signal')
xlabel('sample step t')

subplot(3,1,2)
plot(t,y_hat)
title('Least squares estimated signal')
xlabel('sample step t')

subplot(3,1,3)
hold on
plot(t,y_hat)
plot(t,s)
legend('Estimated signal', 'True signal')
title('True signal overlayed with least squares estimate')
xlabel('sample step t')

print('./homework_2_comparison_2', '-dpng')

