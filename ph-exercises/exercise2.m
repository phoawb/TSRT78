%% 2.3bc

%b)
N = 32;
w1 = 1;
T = 1;
n = 0:(N -1);

y = cos(w1 * n);

[Y, omega] = dtft(y, T, N);
f = omega / (2 * pi);
%f = 0:(N-1);% / N * (1/T);

figure;
plot(f, abs(Y));

%c)

% The leakage comes from the finite representation of the infinite cosine
% signal.

%% 2.5
load("power.mat"); % loads u2 in memory
T = 3.9 * 10 ^ -4;

%a)
y = u2(1:3:end); % every third value of u2

[Y, omega] = dtft(y, 3 * T, length(y));
f = omega / (2 * pi);
[U, omega2] = dtft(u2, T, length(u2));
f2 = omega2 / (2 * pi);

figure;
subplot(121);
plot(omega, abs(Y));
title("DTFT of y");
subplot(122);
plot(omega2, abs(U));
title("DTFT of u2");

figure;
subplot(211);
stem(u2);
title("u2");
subplot(212);
stem(y, "r");
title("y");

%b)
y2 = decimate(u2, 3);

[Y2, omega3] = dtft(y2, 3 * T);

figure;
plot(omega3, abs(Y2));

figure;
subplot(211);
stem(u2);
title("u2");
subplot(212);
stem(y2, "r");
title("low passed y");
