% Define parameters
N = 30;                     % Number of samples
Ts = 1;                     % Sampling interval in seconds
t = (0:N-1)*Ts;             % Time vector
var_e = 0.01;               % Variance of Gaussian noise

% Generate signals
s = sin(t) + sin(1.2*t);    % Signal s(t)
e = sqrt(var_e)*randn(1,N); % Gaussian noise e(t)
y = s + e;                  % Observed signal y(t)

figure;
hold on
tt = 0:0.1:N-1;
plot(tt, sin(tt) + sin(1.2*tt),'b')
stem(t,y,'r')
title('Signal y(t) and sampled points')
xlabel('time')
ylabel('y(t)')
legend('signal', 'N=30 samples')
print('./homework_1_signal', '-dpng')

% Perform the discrete Fourier transform (DFT)
Y = fft(y, N);

% Compute the frequency axis
f = (0:N-1)*(1/(Ts*N));

% Plot the magnitude of the DFT
figure;
stem(f, abs(Y));
title('DFT of y(t)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
print('./homework_1_DFT', '-dpng')

new_N = 2^7;

% Zero-padding to improve frequency resolution
Y_zeropad = fft(y, new_N);

% Compute the new frequency axis for zero-padding
f_zeropad = (0:new_N-1)*(1/(Ts*new_N));

figure;
plot(f_zeropad, abs(Y_zeropad),'b');
title('Zero-padded DFT of y(t)');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
print('./homework_1_DFT_zero', '-dpng')

% DTFT
[Y_DTFT, w] = dtft(y);

figure;
plot(w,abs(Y_DTFT),'b')
xlabel('Angular frequency')
ylabel('|DTFT(y)|')
title('Discrete Time Fourier Transform of y(t)')
print('./homework_1_DTFT', '-dpng')


