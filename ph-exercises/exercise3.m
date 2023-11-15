%%2.25
close all;
%a)

%{
Generate a sinusoid test signal with angular frequency wo = 2pi/5 and
length N = 60. Decimate the signal a factor of two, and compute the
DFT of the resulting signal. Explain the result using the DFT of the original signal.
%}

wo = 2 * pi / 5;
N = 60;
T = 1;
p = 2;
n = (0:N - 1);
x = sin(wo * n);

T2 = T * p;
N2 = N / p;
x2 = x(1:2:N);

X = fft(x, N);
X2 = fft(x2, N2);

ff = (0:(N - 1)) * (1 / (N));
ff2 = (0:(N2 - 1)) * (1 / (T2 * N2));

figure;
subplot(121); stem(n, x);
subplot(122); stem(0:(N2 - 1), x2);

figure;
subplot(211);
stem(ff, abs(X));

subplot(212);
stem(ff2, abs(X2));

%b)

wo = 2 * pi / 5;
N = 60;
T = 1;
p = 3;
n = (0:N - 1);
x = sin(wo * n);

T2 = p * T;
N2 = N / p;
x2 = x(1:p:N);

X = fft(x, N);
X2 = fft(x2, N2);

ff = (0:(N - 1)) * (1 / (N));
ff2 = (0:(N2 - 1)) * (1 / (T2 * N2));

figure;
subplot(121); stem(n, x);
subplot(122); stem(0:(N2 - 1), x2);

figure;
subplot(211);
stem(ff, abs(X));

subplot(212);
stem(ff2, abs(X2));
title("dec 3")

%c)
wo = 2 * pi / 5;
N = 60;
T = 1;
p = 7;
n = (0:N - 1);
x = sin(wo * n);

T2 = T / p;
N2 = round(N / p);
x2 = x(1:p:N);

X = fft(x, N);
X2 = fft(x2, N2);

ff = (0:(N - 1)) * (1 / (N));
ff2 = (0:(N2 - 1)) * (1 / (T2 * N2));

figure;
subplot(121); stem(n, x);
subplot(122); stem(0:(N2 - 1), x2);

figure;
subplot(211);
stem(ff, abs(X));

subplot(212);
stem(ff2, abs(X2));

%% 3.2
close all;
load("sig30.mat"); % load in y

N = length(y);
T = 1;
n = 0:(N-1);

Y = fft(y, N);
ff = (0:(N-1)) * (1/(T * N));

N_zero = 2^16;
Y2 = fft(y, N_zero);
ff2 = (0:(N_zero - 1)) * (1 / (T * N_zero));

figure; 
subplot(211);
plot(ff, abs(Y));
title("FFT of y");
subplot(212);
plot(ff2, abs(Y2));
title("zero-padded fft of y");


%% 4.16

load("sig40.mat") % load s, s1, s2
% w_s1 = 0.1, w_s2 = 0.5 
% s1 = sin(wt), s2 = sin(wt);
% s = s1 + s2;
N = length(s);
n = 0:(N-1);
T = 1; 

[b1, a1] = butter(20, 0.11, 'high'); % Designs the filter with the given order and cutoff, 'type' can be 'low', 'high', 'bandpass', or 'stop'
y1 = filter(b1, a1, s);

[b2, a2] = butter(20, 0.11, "low");
y2 = filter(b2, a2, s);

Y1 = fft(y1);
Y2 = fft(y2);

ff = (0:(N - 1)) * (1 / (T * N));

% PLOT ZONE
figure; 
subplot(511);
plot(n, s1);
title("s1 (0.1 rad/s)");

subplot(512);
plot(n,s2);
title("s2 (0.5 rad/s)");

subplot(513);
plot(n,s);
title("s");


subplot(514);
plot(ff, abs(Y1));

subplot(515);
plot(ff, abs(Y2));

