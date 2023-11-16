close all; clc; clear;

% Sampling time and freq
T_s = 1;
f_s = 1/T_s;

% Number of samples and range
N = 32;
n = 0:1:N-1;

% Angular freq
w_1 = 1;

% Signal
y = cos(n*w_1);

figure;
plot(n,y)

% DTFT
[Y, w] = dtft(y, T_s, N);

figure;
plot(w, abs(Y))

%%
% The reason for the Y looking the way it does is
% because a cut-off in the time domain is a
% convolution with a sinc in the fourier domain
% and the cosine is represented by precisely
% two dirac impulses.
