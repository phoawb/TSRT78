function [X,f]=batchdft(x,S,overlap,Npad)
%Batchwise DFT (spectrogram)
N=length(x);
if nargin<2, S=ceil(N/50); end         % Segment size
if nargin<3, overlap=ceil(0.5*S); end  % Overlapping segment
if nargin<4, Npad=2*S; end             % Zero padding
f=(0:Npad/2-1)/Npad;
w=getwindow(S,'hanning');
M=floor((N-S)/(S-overlap));
for k=1:M;
    Xk=fft([w.*x(k*(S-overlap)+[1:S]);zeros(Npad-S,1)]);
    Xk=Xk(1:Npad/2);
    X(:,k)=1/(w'*w)*real(Xk.*conj(Xk));
end
if nargout==0
    image(f,[1:M]*(S-overlap),256/max(X(:))*X)
    set(gca,'Ydir','normal')
end

