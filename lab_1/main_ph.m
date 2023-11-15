%% Get signal from file and plot data
[y, fSamp] = audioread('data/whistle.wav');

nSamp = size(y, 1);
t = (0:nSamp - 1) / fSamp;

figure;
plot(t, y)
xlabel('time in seconds')
ylabel('recorded signal')

%% Get the energy from the signal

E_tot_time = getEnergy(y, 'time');
E_tot_freq = getEnergy(y, 'frequency');

%compute the energy of the dominant frequency component in the time domain
Fs = 8000; %sampling frequency
N = length(y); %number of samples

% Step 1: Fourier Transform to get frequency spectrum
Y = fft(y);

% Step 2: Find the dominant frequency component
[~, dominantIndex] = max(abs(Y(1:(N / 2 + 1)))); % Only considering the positive frequencies

% Step 3: Isolate the dominant frequency and perform Inverse Fourier Transform
dominantSpectrum = zeros(size(Y)); % Initialize with zeros
dominantSpectrum(dominantIndex) = Y(dominantIndex); % Retain only the dominant frequency
dominantSpectrum(end - dominantIndex + 2) = Y(end - dominantIndex + 2); % Retain the symmetric component
dominantSignal = real(ifft(dominantSpectrum)); % Inverse FFT

% Step 4: Calculate the energy of the dominant frequency component
E_dom_time = sum(dominantSignal .^ 2)
