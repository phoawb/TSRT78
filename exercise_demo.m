%% Basic parameters
N = 1000;
n = 1:N-1;
sigma_noise = 0.04;
w = 0.2;
y = sin(w*n) + sigma_noise*randn(1,N-1);
figure, plot(n,y);

%% AR(2)
ar_mod = ar(y,2);

%% Periodogram
non_parametric_mod = etfe(y);
figure, bode(non_parametric_mod)
hold on, bode(spa(y', 200));

%% Pole zero diagram
figure, zplane(1, ar_mod.A)

% get the argument:
wT = angle(poles);

T = 1; % in this ex
w_est = wT/T;

% Get residuals pe gives us prediction error
eps = pe(ar_mod, y');
figure;
subplot(1,2,1), scatter(1:length(eps), eps), ylabel('\epsilon')
subplot(1,2,2), qqplot(eps) % see if the error is normally distrib (if match red line)

% calculate loss V_N(theta)
loss = 1/length(eps) * sum(eps.^2);

Ree = covf(eps, 20);
figure, plot(0:20-1, Ree/Ree(1))
ylabel('R_{ee}(\tau)')
xlabel('\tau')