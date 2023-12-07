
load("lunarmodule.mat")

A = [1 1; 0 1];
B = [1; 1];
C = [1 0];
Qtil = 1;
R = 500;

[yhat, xhat] = kalmanfilt(y, A, B, C, Qtil, R);
disp(yhat);
disp(xhat);