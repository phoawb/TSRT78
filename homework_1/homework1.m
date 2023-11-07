% Time instances
t = 1:0.1:100;
% Sampling interval
T_s = 1;
N = 30;

s = sin(t) %+ sin(1.2 * t);

figure;
plot(t, s);

disp(t);
