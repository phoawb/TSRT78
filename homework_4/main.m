clear; close all; clc;
load Aug11.mat

meas = ykf; % get measurement

%ykf=detrend(ykf);
n_meas = length(meas); % sequence length
time = 1:300; % time steps

figure;
plot(1:n_meas, meas, 'linewidth', 2); % plot measurement
title("measured signal y")

% dim of the parameter vector
dim = 2;

[tht, Pt] = special_kf(ykf, dim);

phi = ones(2, 300);
phi(2, :) = 1:300; % initialize the regression vector

s_est = diag(phi' * tht)'; % get the estimates of s


figure;
plot(1:n_meas, tht);
title("theta")
legend("a_t", "b_t")

figure;
plot(1:n_meas, s_est);
title("estimated signal s")

figure;
plot(1:n_meas, meas, 'linewidth', 2); hold on; % plot measurement
plot(1:n_meas, s_est, 'r', 'linewidth', 2); % plot estimated signal
legend('Measurements', 'Estimated Signal');
title('Signal Estimation using adaptive Kalman Filter');
xlabel('Time');
ylabel('Signal Value');


function [tht,Pt]=special_kf(y,na,Q)
%RLS for adaptive AR parameter estimation
if nargin<3; Q=(1e-5)*eye(na); end

% init variables 
R = 1;
N=length(y);
th=zeros(na,1);
epsi = zeros(1, N);
tht = zeros(na, N);
P=100*eye(na);
Pt = zeros(na, na, N);

for t=1:N
    phi = [1 t]';
    K= P * phi ./ (R + phi' * P * phi);

    %update & store
    P = P - K * phi' * P + Q;
    epsi(:,t) = y(t) - phi' * th;
    th = th + K * epsi(:,t);
    tht(:,t) = th; 
    Pt(:,:,t) = P;
end
end


 




