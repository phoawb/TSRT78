% Define parameters
N = 30;                     % Number of samples
Ts = 1;                     % Sampling interval in seconds
t = (0:N-1)*Ts;             % Time vector
var_e = 0.01;               % Variance of Gaussian noise

% Generate signals
s = sin(t) + sin(1.2*t);    % Signal s(t)
e = sqrt(var_e)*randn(1,N); % Gaussian noise e(t)
y = s + e;                  % Observed signal y(t)


[Y_dtft, omega] = dtft(y, T, N);

figure;
plot(omega, abs(Y_dtft));
title("Magnitude of the DTFT of y(t)")
ylabel("|DTFT{y(t)}|")
xlabel("omega (rad/sec)")

% Perform the discrete Fourier transform (DFT)
Y = fft(y);

% Compute the frequency axis
f = (0:N-1)*(1/(Ts*N));

% Plot the magnitude of the DFT
figure;
stem(f, abs(Y));
title('Magnitude of the DFT of y(t)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');

N_zeropad = 2^10;

% Zero-padding to improve frequency resolution
Y_zeropad = fft(y, N_zeropad);  % Zero-pad to the next power of 2 greater than N

% Compute the new frequency axis for zero-padding
f_zeropad = (0:N_zeropad-1)*(1/(Ts*N_zeropad));

% Plot the magnitude of the zero-padded DFT
figure;
plot(f_zeropad, abs(Y_zeropad));
title('Zero-padded DFT of y(t) (N = 1024)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
