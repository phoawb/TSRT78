%% 6.5
close all; clear;
%a)

N = 300;
n = 0:(N - 1);
w = 0.1;

s = sin(w * n);

figure;
plot(n, s);

t1 = ar(s, 2);
%t1.a gives parameters
r1 = roots(t1.a);
fr1 = angle(r1);

figure; bode(t1);

figure; zplane(1, t1.a)

%b)
sigma_noise = 0.01
s = s + sigma_noise * randn(1, N);
figure;
plot(n, s);

t2 = ar(s, 2);
r2 = roots(t2.a);
fr2 = angle(r2)

figure; bode(t2)

%% 6.18

%a)
load("sig60.mat"); % loads y into memory
y = y';
N = length(y);
N_est = round((2 / 3) * N);
N_validate = round((1 / 3) * N);

y_est = y(1:N_est);
y_validate = y(N_est+1:end);

t1 = ar(y_est, 1);
t2 = ar(y_est, 2);
t3 = ar(y_est, 3);

le(1) = (1/N_est) * sum(pe(t1, y_est') .^2);
le(2) = (1/N_est) * sum(pe(t2, y_est') .^2);
le(3) = (1/N_est) * sum(pe(t3, y_est') .^2);

t1 = ar(y_validate, 1);
t2 = ar(y_validate, 2);
t3 = ar(y_validate, 3);

lv(1) = (1/N_validate) * sum(pe(t1, y_validate') .^2);
lv(2) = (1/N_validate) * sum(pe(t2, y_validate') .^2);
lv(3) = (1/N_validate) * sum(pe(t3, y_validate') .^2);

figure;
plot(1:3, le, "-")
hold on
plot(1:3, lv, "--")

