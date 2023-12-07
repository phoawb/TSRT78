close all; clc; clear;
[vowel_o_raw, fs] = audioread('data/ooo.wav');

% ---------- PRE-PROCESSING -------------------------
vowel_o_raw = detrend(vowel_o_raw);
N = size(vowel_o_raw,1); % number of samples
t = (0:N-1)/fs; % time vector in seconds
figure;clf;
plot(t, vowel_o_raw)
title('ooo raw')
xlabel('time in seconds')
ylabel('recorded signal') % axis description is important!

vowel_o = vowel_o_raw(floor(8000*6.5)+1:floor(8000*8.5)); % cut out best 2 seconds

N = size(vowel_o,1); % number of samples
t = (0:N-1)/fs; % time vector in seconds

figure;clf;
plot(t, vowel_o)
axis([])
xlabel('time in seconds')
ylabel('recorded signal') % axis description is important!

% ---------- EST MODEL ORDER FROM GRAPH -------------
% 6*2 peaks = order 12 appropriate for AR
VOWEL_O = fft(vowel_o,N);
ff = (0:(N-1))*(fs/N);
figure; clf;
plot(ff,abs(VOWEL_O));      
title("Fourier transform of the a sound")
xlabel("frequency (Hz)")
ylabel("amplitude of Fourier transform")

% ---------- PARTITION DATA -------------------------
% 2/3 for estimation 1/3 for validation
vowel_o_est = vowel_o(1:2*floor(N/3));      % Estimated data
vowel_o_val = vowel_o(2*floor(N/3)+1:N);    % Validation data

% ---------- PRE-defined LOSS function plot ---------
nmax = 20;
[n, W_cv] = arordercv(vowel_o_est, vowel_o_val, nmax);
disp('best order')
disp(n)

%% ---------- HOME-MADE SOLUTION ---------------------
[W, Uaic, Ubic] = arorder(vowel_o, nmax);

figure;
plot(3:nmax, W_cv(3:end), "-o");
title("Cross validation loss function")
xlabel("AR model order")
ylabel("loss")
xticks(3:nmax)

figure;
subplot(2,1,1);
plot(1:nmax,W,'-o', 1:nmax,Uaic, '-.dr', 1:nmax, Ubic, '--xb')
title('Loss functions as a function of AR(n)')
legend('Loss function', 'AIC', 'BIC')
xlabel('AR model order')
ylabel("loss")
xticks(1:nmax)
subplot(2,1,2);
plot(1:nmax,W,'-o', 1:nmax,Uaic, '-.dr', 1:nmax, Ubic, '--xb')
xlim([3 20])
title('Loss functions as a function of AR(n)')
legend('Loss function', 'AIC', 'BIC')
xlabel('AR model order')
ylabel("loss")
xticks(1:nmax)


%% ---------- SPECTRUM COMPARISON ------------------
data_est = iddata(vowel_o_est, [], 1/fs);
data_val = iddata(vowel_o_val, [], 1/fs);
vowel_o_spect = etfe(data_val, 200);
ar_3 = ar(data_est,3);
ar_16 = ar(data_est,16);
ar_30 = ar(data_est,30);

figure;clf;
bode(ar_3, 'g-', ar_16, 'c-', ar_30, 'b-', vowel_o_spect, 'r');
legend('ar(3)', 'ar(16)', 'ar(30)', 'periodogram')

%% ---------- COMPARE and PREDICT ------------------
figure;
subplot(3,1,1)
YP = predict(ar_3, data_val);
compare(data_val, YP);
subplot(3,1,2)
YP = predict(ar_16, data_val);
compare(data_val, YP);
subplot(3,1,3)
YP = predict(ar_30, data_val);
compare(data_val, YP);

%% ------------------ Simulation of model ------------------
mo18 = ar(vowel_o,3); 
e_vec = filter(mo18.a,1,vowel_o); % E = (a_1 + ...) * Y

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
plot(1:length(vowel_o),vowel_o);
hold on;
plot(1:length(vowel_a_out),vowel_a_out);
hold off


