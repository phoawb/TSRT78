close all
% Generate a 10^4 long white noise sample sequence
N = 10^4;
u = randn(1,N);

tt = 0:N-1;

% Channel from impulse response
% h(t) = 0.9δ(t−3) − 0.4δ(t−4) + 0.2δ(t−7)
h = zeros(1,8);
h(4) = 0.9;
h(5) = -0.4;
h(8) = 0.2;

% filtered_signal = conv(u,h,"same");
y = filter(h,1,u);

% Apply the LMS algorithm
nb = length(h);        % Model order
nk = 2;                % Adjust as needed
mu = 0.015;            % step-length
lambda = 0.005;        % leakage factor

[th, ~, ~] = MyLMS_2(y', u', nb, nk, mu, lambda);


figure;
subplot(3,1,1)
plot(tt,u)

subplot(3,1,2)
stem(1:8, h)

subplot(3,1,3)
plot(tt,y)


figure;
% Plotting the evolution of filter coefficients
plot(th)
title('Evolution of LMS Estimated Filter Coefficients')
legend('Coefficient 1', 'Coefficient 2', 'Coefficient 3', 'Coefficient 4', 'Coefficient 5', 'Coefficient 6', 'Coefficient 7', 'Coefficient 8')
xlabel('Number of Samples')
ylabel('Coefficient Value')

% Assuming 'th' contains the parameter estimates from the LMS algorithm
final_estimate = th(end, :);

% Actual parameters that generated the noise sequences
actual_params = [0, 0, 0, 0.9, -0.4, 0, 0, 0.2];

% Create the plot
figure;
hold on;
stem(1:length(final_estimate), final_estimate, 'o', 'MarkerFaceColor', 'blue', 'DisplayName', 'LMS-param.');
stem(1:length(actual_params), actual_params, '*', 'MarkerFaceColor', 'red', 'DisplayName', 'Param. for the noise sequence');
hold off;

% Add details
legend('show');
xlabel('Parameter');
ylabel('Value');
title('Parameters at the Last Iteration and Parameters that Generate the Noise Sequences');
grid on; % Add a grid for better readability

%{
QUESTION 1:
* Geometry & distance: 
    Yes, by measuring the distance and speed of sound.
    Drawback is that we don't consider the inherent delay of
    the electronics.
* Visual examinatinon of the signals by using Matlab plot
    YES, the delay is visible because
    both samples just start from 0.
* Computing with xcorr(y,u,M)
    Yes it is possible, however expensive
* Estimating channel parameter
    It is possible, by inspecting where the coeffiecients
    begin. Assuming nk is zero. Otherwise subtract nk.

QUESTION 2:
* I the step length is too large, it will never converge.
    Takes longer time, more expensive for smaller step len.
    More precise, though.
`
QUESTION 3:
* The resulting error, |yhat - y|

QUESTION 4:
* see file

%}

%%
[A_rec_white,B_rec_white]=play_and_rec_noise("/dev/cu.usbmodem11101", 'white')
save white.mat A_rec_white B_rec_white 
%%
[A_rec_sine,B_rec_sine]=play_and_rec_noise("/dev/cu.usbmodem11101", 'sine')
save sine.mat A_rec_sine B_rec_sine 
%%
[A_rec_multisine,B_rec_multisine]=play_and_rec_noise("/dev/cu.usbmodem11101", 'multisine')
save multisine.mat A_rec_multisine B_rec_multisine 
%%
[A_rec_chirp,B_rec_chirp]=play_and_rec_noise("/dev/cu.usbmodem11101", 'chirp')
save chirp.mat A_rec_chirp B_rec_chirp 

%% Load recorded signals from microphones A and B
chirp_a = load('chirp.mat').A_rec_chirp;
chirp_b = load('chirp.mat').B_rec_chirp;
sine_a = load('sine.mat').A_rec_sine;
sine_b = load('sine.mat').B_rec_sine;
multisine_a = load('multisine.mat').A_rec_multisine;
multisine_b = load('multisine.mat').B_rec_multisine;
white_a = load('white.mat').A_rec_white;
white_b = load('white.mat').B_rec_white;

% Normalize signals and detrend
chirp_a = chirp_a/max(chirp_a);
chirp_a = detrend(chirp_a);

chirp_b = chirp_b/max(chirp_b);
chirp_b = detrend(chirp_b);

sine_a = sine_a/max(sine_a);
sine_a = detrend(sine_a);

sine_b = sine_b/max(sine_b);
sine_b = detrend(sine_b);

multisine_a = multisine_a/max(multisine_a);
multisine_a = detrend(multisine_a);

multisine_b = multisine_b/max(multisine_b);
multisine_b = detrend(multisine_b);

white_a = white_a/max(white_a);
white_a = detrend(white_a);

white_b = white_b/max(white_b);
white_b = detrend(white_b);
%% Define sampling frequency and time

fs = 8000; % [Hz]
Ts = 1/fs;

%% nk With xcorr
%{
we consistently get 12, only differing at sine which gave 10'
%}
close all

maxLag = 70;
[R,lags] = xcorr(chirp_a,chirp_b,maxLag);
figure;
subplot(4,1,1);
plot(lags,R); xlabel('Lags'); title('Chirp correlation for different lags'); 

[R,lags] = xcorr(sine_a,sine_b,maxLag);
subplot(4,1,2);
plot(lags,R); xlabel('Lags'); title('Sine correlation for different lags'); 

[R,lags] = xcorr(multisine_a,multisine_b,maxLag);
subplot(4,1,3);
plot(lags,R); xlabel('Lags'); title('Multisine correlation for different lags'); 

[R,lags] = xcorr(white_a,white_b,maxLag);
subplot(4,1,4);
plot(lags,R); xlabel('Lags'); title('White correlation for different lags'); 

%% nk with high model order
%{
Answer, we got 11 on multisine and white, but the other
ones were hard to gauge, so nk is 10
%}

nk = 2;                         % time-delay in sample
nb = 100;                       % nb + nk = model order 
mu = 0.015;                     % step-length
lambda = 0.005;                 % leakage factor

[th_chirp,sHat_chirp,err_chirp] = MyLMS_2(chirp_a, chirp_b, nb, nk, mu, lambda);
[th_sine,sHat_sine,err_sine] = MyLMS_2(sine_a, sine_b, nb, nk, mu, lambda);
[th_multisine,sHat_multisine,err_multisine] = MyLMS_2(multisine_a, multisine_b, nb, nk, mu, lambda);
[th_white,sHat_white,err_white] = MyLMS_2(white_a, white_b, nb, nk, mu, lambda);

figure;
subplot(4,1,1)
stem(1:length(th_chirp(end,:)), th_chirp(end,:), 'o', 'MarkerFaceColor', 'blue', 'DisplayName', 'LMS-param.');
subplot(4,1,2)
stem(1:length(th_sine(end,:)), th_sine(end,:), 'o', 'MarkerFaceColor', 'blue', 'DisplayName', 'LMS-param.');
subplot(4,1,3)
stem(1:length(th_multisine(end,:)), th_multisine(end,:), 'o', 'MarkerFaceColor', 'blue', 'DisplayName', 'LMS-param.');
subplot(4,1,4)
stem(1:length(th_white(end,:)), th_white(end,:), 'o', 'MarkerFaceColor', 'blue', 'DisplayName', 'LMS-param.');

%% ========================= Determine best nb ============================

%% ---------------- With AIC/BIC/loss plots
figure(8);
maxOrder = 100;
subplot(4,1,1);
arorder(chirp_a,maxOrder);
legend('Minimal loss function',...
       'Akaikes information criterion (AIC)',...
       'Akikes information criterion B (BIC)');
title('Chirp model order');

subplot(4,1,2);
arorder(sine_a,maxOrder);
legend('Minimal loss function',...
       'Akaikes information criterion (AIC)',...
       'Akikes information criterion B (BIC)');
title('Sine model order');

subplot(4,1,3);
arorder(multisine_a,maxOrder);
legend('Minimal loss function',...
       'Akaikes information criterion (AIC)',...
       'Akikes information criterion B (BIC)');
title('Multisine model order');

subplot(4,1,4);
arorder(white_a,maxOrder);
legend('Minimal loss function',...
       'Akaikes information criterion (AIC)',...
       'Akikes information criterion B (BIC)');
title('White model order');

%{
ANSWER:
Looking at the loss function plots we
observe:

multisine order: 40 (has many peaks)
chirp order: 6 (peaks hard to count)
sine order: 2 (sine only has two peaks)
white order: 43 (has infinite peaks)
%}
%% Number of peaks
fs = 8000;

N = length(chirp_b)*10;
ff = (0:(N-1)) * (fs/N);

figure;
subplot(2,1,1)
plot(ff, abs(fft(chirp_b,N)))

subplot(2,1,2)
plot(ff, abs(fft(multisine_b,N)))

% very many peaks on chirp spectra, impossible to count

%% LMS expectation

%% Estimate model parameters with MyLMS

nk = 11;                        % time-delay in sample
nb = 43;                        % nb + nk = model order 
mu = 0.015;                     % step-length
lambda = 0.005;                 % leakage factor

[th_chirp,sHat_chirp,err_chirp] = MyLMS_2(chirp_a, chirp_b, nb, nk, mu, lambda);
[th_sine,sHat_sine,err_sine] = MyLMS_2(sine_a, sine_b, nb, nk, mu, lambda);
[th_multisine,sHat_multisine,err_multisine] = MyLMS_2(multisine_a, multisine_b, nb, nk, mu, lambda);
[th_white,sHat_white,err_white] = MyLMS_2(white_a, white_b, nb, nk, mu, lambda);

figure(1);
subplot(4,1,1);
plot(th_chirp); xlabel('Iteration'); title('Chirp model parameters');
subplot(4,1,2);
plot(th_sine); xlabel('Iteration'); title('Sine model parameters');
subplot(4,1,3);
plot(th_multisine); xlabel('Iteration'); title('Multisine model parameters');
subplot(4,1,4);
plot(th_white); xlabel('Iteration'); title('White model parameters');

%% ======================= 3. Determine energy of sHat ====================

chirp_energy = sum(abs(chirp_a).^2);
sine_energy = sum(abs(sine_a).^2);
multisine_energy = sum(abs(multisine_a).^2);
white_energy = sum(abs(white_a).^2);

sHat_chirp_energy = sum(abs(sHat_chirp).^2);
sHat_sine_energy = sum(abs(sHat_sine).^2);
sHat_multisine_energy = sum(abs(sHat_multisine).^2);
sHat_white_energy = sum(abs(sHat_white).^2);

sHat_chirp_energy/chirp_energy             % 0.0071
sHat_sine_energy/sine_energy               % 0.0037
sHat_multisine_energy/multisine_energy     % 0.5045
sHat_white_energy/white_energy             % 0.6358

%{
Best estimates for sine and chirp
not so good for multisine and white
%}

%%  Compare real signal and its estimation
close all;

yHat_chirp = chirp_a - sHat_chirp;
yHat_sine = sine_a - sHat_sine;
yHat_multisine = multisine_a - sHat_multisine;
yHat_white = white_a - sHat_white;

clf;
subplot(4,1,1);
plot(chirp_a); hold on; plot(yHat_chirp); title('Chirp: real and estimation');
legend('y','yHat');
subplot(4,1,2);
plot(sine_a); hold on; plot(yHat_sine); title('Sine: real and estimation');
legend('y','yHat');
subplot(4,1,3);
plot(multisine_a); hold on; plot(yHat_multisine); title('Multisine: real and estimation');
legend('y','yHat');
subplot(4,1,4);
plot(white_a); hold on; plot(yHat_white); title('White: real and estimation');
legend('y','yHat');

%{
QUESTION 4:

Sine and chirp are very close, can also guessed from the
energy ratio, which was low.

White is pretty baddly estimated, which was also expected.

Multisine was also difficult to dynamically model, pretty
bad actually. Only slightly better than white noise.

QUESTION 5:

It was easiest to supress noise from the sine and chirp.
The sine was easy to model since it is a very simple signal with
few frequencies, the chirp was also very easy to model, for some
reason.
%}