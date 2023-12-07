close all; clear; clc;

% sinusoid coefficients
A = 1 / sqrt(2);
xi = 0.1;
phi = pi / 2;

N = 1000; % number of samples
var = 0.25; % variance of noise
e = sqrt(var) * randn(1, N); % zero-mean Gaussian noise

tt = 1:N;
s = A * sin(xi * tt + phi); %sinusoid
y = s + e; % noisy measurement

% sin(a + b) = sin(a)cos(b) + cos(a)sin(b)
% sin(xi * tt) & cos(xi * tt) are known so that becomes our data matrix
% b = [Chat; Dhat] where C = A * cos(phi) and D = A * sin(phi)
% (X' * X) * b = X' * y

% Constructing the design matrix X
X = [sin(xi * tt); cos(xi * tt)]';

% Least Squares Estimation
b = (X' * X) \ X' * y';

% Extract estimates
Chat = b(1);
Dhat = b(2);

% Estimate amplitude and phase shift
Ahat = sqrt(Chat ^ 2 + Dhat ^ 2);
phihat = atan2(Dhat, Chat);

% Display results
disp(['Estimated Amplitude: ', num2str(Ahat)]);
disp(['Estimated Phase Shift: ', num2str(phihat)]);

% PLOT ZONE
figure;
subplot(211); plot(tt, s);
title("Sinusoid")
subplot(212); plot(tt, y);
title("Noisy measurements");

figure;
subplot(311); plot(tt, s);
title("Original Sinusoid")
subplot(312); plot(tt, y);
title("Noisy measurements");
subplot(313); plot(tt, Ahat * sin(xi * tt + phihat));
title("Estimated Sinusoid");

figure;
plot(tt, s, 'b', tt, Ahat * sin(xi * tt + phihat), '--r');
title("Original Sinusoid vs Estimated Sinusoid");
legend("Original Sinusoid", "Estimated Sinusoid");
