close all; clc; clear;
[vowel_a_raw, fs] = audioread('data/aaa.wav');

% ---------- PRE-PROCESSING -------------------------
vowel_a_raw = detrend(vowel_a_raw);
N = size(vowel_a_raw,1); % number of samples
t = (0:N-1)/fs; % time vector in seconds
figure;clf;
plot(t, vowel_a_raw)
title('aaa raw')
xlabel('time in seconds')
ylabel('recorded signal') % axis description is important!

vowel_a = vowel_a_raw(floor(8000*1.7)+1:floor(8000*3.7)); % cut out best 2 seconds

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
disp('best order')
disp(n)

%% ---------- HOME-MADE SOLUTION ---------------------
[W, Uaic, Ubic] = arorder(vowel_a, nmax);

figure;
plot(3:nmax, W_cv(3:end), "-o");
title("Cross validation loss function")
xticks(3:nmax)

figure;
subplot(2,1,1);
plot(1:nmax,W,'-o', 1:nmax,Uaic, '-.dr', 1:nmax, Ubic, '--xb')
title('Loss functions as a function of AR(n)')
legend('Loss function', 'AIC', 'BIC')
xlabel('n')
xticks(1:nmax)
subplot(2,1,2);
plot(1:nmax,W,'-o', 1:nmax,Uaic, '-.dr', 1:nmax, Ubic, '--xb')
xlim([3 16])
title('Loss functions as a function of AR(n)')
legend('Loss function', 'AIC', 'BIC')
xlabel('n')
xticks(1:nmax)


%% ---------- SPECTRUM COMPARISON ------------------
data_est = iddata(vowel_a_est, [], 1/fs);
data_val = iddata(vowel_a_val, [], 1/fs);
vowel_a_spect = etfe(data_val, 200);
ar_3 = ar(data_est,3);
ar_18 = ar(data_est,18);
ar_30 = ar(data_est,30);

figure;clf;
bode(ar_3, 'g-', ar_18, 'c-', ar_30, 'b-', vowel_a_spect, 'r');
legend('ar(3)', 'ar(18)', 'ar(30)', 'periodogram')

%% ---------- COMPARE and PREDICT ------------------
figure;
subplot(3,1,1)
YP = predict(ar_3, data_val);
compare(data_val, YP);
subplot(3,1,2)
YP = predict(ar_18, data_val);
compare(data_val, YP);
subplot(3,1,3)
YP = predict(ar_30, data_val);
compare(data_val, YP);

%% ------------------ Simulation of model ------------------
mo18 = ar(vowel_a,18); % 18 because we have 9 clear peaks
e_vec = filter(mo18.a,1,vowel_a); % E = (a_1 + ...) * Y

f_min = 20;
f_max = 100;

r = covf(e_vec,f_max);
[r_max,r_i] = max(r(f_min:end));
f_period = f_min + r_i;

pulsetrain = zeros(1, N);
pulsetrain(1:f_period:end)=sqrt(r_max); % speech synthesis hack

vowel_a_out = filter(1, mo18.a, pulsetrain); % Y = 1/(a_1 + ...) * E
sound(50*vowel_a_out, fs);

figure; clf;
plot(1:length(vowel_a),vowel_a);
hold on;
plot(1:length(vowel_a_out),vowel_a_out);
hold off


