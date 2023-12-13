function [Phi,f]=sig2periodogram(x,T)
if nargin<2, T=1; end
N=length(x);
X=fft(x);
Phi=X.*conj(X);
Phi=Phi(1:round(N/2)+1)/N*T;
f=(0:round(N/2))/N/T;