close all; clc; clear;
[vowel_a_raw, fs] = audioread('Sound/A2.wav'); 

% ---------- PRE-PROCESSING -------------------------
vowel_a_raw = detrend(vowel_a_raw);
vowel_a = vowel_a_raw(floor(8000*3.27)+1:floor(8000*5.27)); % cut out best 2 seconds

N = size(vowel_a,1); % number of samples
t = (0:N-1)/fs; % time vector in seconds

figure;clf;
plot(t, vowel_a)
axis([])
xlabel('time in seconds')
ylabel('recorded signal') % axis description is important!

% ---------- EST MODEL ORDER FROM GRAPH -------------
% 6*2 peaks = order 12 appropriate for AR
VOWEL_A = fft(vowel_a,N);
ff = (0:(N-1))*(fs/N);
figure; clf;
plot(ff,abs(VOWEL_A));      

% ---------- PARTITION DATA -------------------------
% 2/3 for estimation 1/3 for validation
vowel_a_est = vowel_a(1:2*floor(N/3));      % Estimated data
vowel_a_val = vowel_a(2*floor(N/3)+1:N);    % Validation data

% ---------- PRE-defined LOSS function plot ---------
nmax = 20;
[n, W_cv] = arordercv(vowel_a_est, vowel_a_val, nmax);
disp("BEST ORDER = ")
disp(n)

% ---------- HOME-MADE SOLUTION ---------------------
lams = [];
Wv_cv = zeros(1, nmax);

[W, Uaic, Ubic] = arorder(vowel_a, nmax);


figure;
plot(3:nmax, W_cv(3:end), "-o");
title("Cross validation loss function")
xticks(3:nmax)

figure;
plot(1:nmax,W,'-o', 1:nmax,Uaic, '-.dr', 1:nmax, Ubic, '--xb')
title('Loss functions as a function of AR(n)')
legend('Loss function', 'AIC', 'BIC')
xlabel('n')
xticks(1:nmax)

% ---------- SPECTRUM COMPARISON of chosen order --
data_est = iddata(vowel_a_est, [], 1/fs);
data_val = iddata(vowel_a_val, [], 1/fs);
vowel_a_spect = etfe(data_val, 200);
order_n = 20;

figure;clf;
bode(ar(data_est,order_n), 'b', vowel_a_spect, 'r');

% ---------- COMPARE and PREDICT ------------------
AR = ar(data_est,order_n);
predict(AR, vowel_a_val);
figure;
compare(data_val, AR);


% Assuming AR = ar(y_est, 9) is your model
YP = predict(AR, vowel_a_val);

% Plot the results
t_val = (1:length(vowel_a_val)); % Time or sample indices for validation data
plot(t_val, vowel_a_val, 'b', t_val, YP, 'r--');
legend('Actual Validation Data', 'Predicted Data');
xlabel('Time or Samples');
ylabel('Signal Amplitude');
title('Comparison of Validation Data and AR Model Predictions');


