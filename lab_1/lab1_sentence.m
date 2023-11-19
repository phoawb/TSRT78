close all; clc; clear;
fn = 4000;
fs = 8000;
% ---- Speech encoding as in GSM ---------------------

[sen, fSamp] = audioread('data/sen2.wav');  % extract data
sen = [sen; zeros(1,16000-length(sen))'];   % zero-pad end with 0's
%sound(sen,fSamp);
N = size(sen,1);                        % number of samples
t = (0:N-1)/fSamp;                      % time vector in seconds

%% Split up sentence with corresponding AR models in 100 chunks
order = 15;
sen_div=zeros(100,160);
for i=1:100
   sen_div(i,:)=sen(160*(i-1)+1:160*(i));
end

sen_ar = zeros(100, order+1);
for i=1:100
   % Get AR model
   temp = ar(sen_div(i,:), order);

   % Check for pole stability (above order 33 is unstable sometimes)
   [R,P,K] = residue(1,temp.a);
   for j=1:length(P)
       dist = norm(P(j));
       if dist > 1
           P(j) = 1/P(j);
       end  
   end 
   [~, new_a] = residue(R,P,K);

   % Set the ar model in list
   sen_ar(i,:) = new_a; 
end

%% Get the covariance function for the middle chunk i=50
chunk_n = 50;
e_vec = filter(sen_ar(chunk_n,:),1,sen_div(chunk_n,:)); %
r = covf(e_vec',100);
figure(7)
plot(r);title('Covariance Function');xlabel('t');
print(7,'convf.eps','-depsc','-loose');
    
%% 
sen_div_hat = zeros(100,160);
f_min = 20;
f_max = 100;

for i=1:100
    e_vec = filter(sen_ar(i,:),1,sen_div(i,:)); 
    r = covf(e_vec', f_max);              % Covariance matrix 
    [r_max,r_index] = max(r(f_min:end));  % skip first f_min
    r_index = r_index + f_min;            % adjust for offset

    pulsetrain = zeros(1,160);               % declare vector
    pulsetrain(1:r_index:end)=sqrt(r_max);   % speech synthesis hack
    sen_div_hat(i,:) = filter(1,sen_ar(i,:),pulsetrain);
end

% Detrend the estimated signals
sen_div_hat_detrend = zeros(1, N);
for i=1:100
    sen_div_hat_detrend(160*(i-1)+1:160*(i)) = detrend(sen_div_hat(i,:));
end

%% Figure zone

tt = (0:(N-1))/fs;

SEN = fft(sen);
SEN_DIV_HAT_DETREND = fft(sen_div_hat_detrend);

figure(4);
subplot(2,1,1)
plot(tt, sen_div_hat_detrend);
title('Parametric Sentence');
xlabel('time (s)');
ylabel('Amplitude')

subplot(2,1,2)
plot(tt, sen);
title('Non-Parametric Sentence');
xlabel('time (s)');
ylabel('Amplitude')
print(4,'sent_t.eps','-depsc','-loose');

ff = (0:N-1)*(fs/N);

figure(5); 
subplot(2,1,1);
plot(ff, abs(SEN_DIV_HAT_DETREND));
title('Parametric Sentence');
xlabel('frequency (Hz)');
ylabel('Amplitude')
subplot(2,1,2)
plot(ff, abs(SEN));
title('Non Parametric Sentence');
xlabel('frequency (Hz)');
ylabel('Amplitude')
print(5,'sen_f.eps','-depsc','-loose');
%sound(y,8000)
%% Play the parametric recording

sound(sen_div_hat_detrend,fs);

%%

function [ new_a ] = pole_stable(a)
% Make sure are the poles are within the
% unit circle

end
