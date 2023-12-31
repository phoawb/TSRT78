function [yk,mpred]=lsspredict(m,y,kstep)
%Prediction for the innovation model
if any(m.Q~=m.R) | any(m.Q~=m.S), 
    error('Model not on innovation form'), 
end
yk=zeros(size(y)+[kstep 0]);
tmp=m.C*m.A^kstep;
xhat=zeros(size(m.A,1),1); 
for l=1:length(y)
   yk(l+kstep)=tmp*xhat;
   xhat=(m.A-m.B*m.C)*xhat+m.B*y(l);
end
% State space model for predictor
mpred.A=m.A-m.B*m.C;  
mpred.B=m.B; mpred.C=tmp; mpred.D=0;

