close all; clear; clc;
load power

figure;
subplot(2,1,1);
hold on;
stem(u2, 'b')
title('base u2')

subplot(2,1,2);
y = u2(1:3:end);
stem(y, 'r')
title('sampled u2')

figure;
[U, w1] = dtft(u2);
[Y, w2] = dtft(y);
hold on;
title('Comparison of DTFT')
plot(w1, abs(U), 'b');
plot(w2, abs(Y), 'r');
legend('U','Y');

%% low sampled u2
figure;
subplot(2,1,1);
hold on;
stem(u2, 'b')
title('u2')

subplot(2,1,2);
y = decimate(u2, 3);
stem(y, 'r')
title('decimated u2')

figure;
[U, w1] = dtft(u2);
[Y, w2] = dtft(y);
hold on;
plot(w1, abs(U), 'b');
plot(w2, abs(Y), 'r');
legend('U','Y');

