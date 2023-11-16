[aaa, fSamp] = audioread('Sound/A2.wav'); 
%sound(aaa,fSamp); 
nSamp = size(aaa,1);
t = (0:nSamp-1)/fSamp; 

figure(1);clf();
plot(t, aaa)
xlabel('time in seconds')
ylabel('recorded signal') % axis description is important!

aa = aaa(floor(8000*3.27)+1:floor(8000*5.27));

nSamp = size(aa,1); % number of samples
t = (0:nSamp-1)/fSamp; % time vector in seconds

figure(2);clf();
plot(t, aa)
xlabel('time in seconds')
ylabel('recorded signal') % axis description is important!

figure(3);clf();
plot(t, aa)
axis([])
xlabel('time in seconds')
ylabel('recorded signal') % axis description is important!

%% Count number of peaks
AA = fft(aa);

%% Cross validation

T= 1/8000;
N=length(aa);

edata = iddata(aa(1:2*floor(N/3)),[],T); % Estimated data
vdata = iddata(aa(2*floor(N/3)+1:N),[],T); % Validation data

for n=1:20
end