%% Get signal from file and plot data
[y, fSamp] = audioread('data/whistle.wav');
    
nSamp = size(y, 1);
t = (0:nSamp - 1) / fSamp;

figure;
plot(t, y)
xlabel('time in seconds')
ylabel('recorded signal')

%% Get the energy from the signal

E_t = getEnergy(y, 'time');
