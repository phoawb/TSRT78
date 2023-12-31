function [h,k]=lss2pulseresponse(m,kmax)
%State space model to pulse response h(k) conversion
k=0:kmax;
h=zeros(size(k));
h(1)=m.D; 
tmp=m.B;
for l=1:kmax
   h(l+1)=m.C*tmp;
   tmp=m.A*tmp;
end
if nargout==0, plot(k,h), end