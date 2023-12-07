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
playsound("sine")

%%
play_and_rec_noise("/dev/cu.usbmodem11101", 'white')
