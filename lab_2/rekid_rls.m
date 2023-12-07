function [th, yhat] = rekid_rls(y,u,N,lambda)

ny=length(y);
th=zeros(N+1,ny+1); yhat=zeros(ny,1);
P=eye(N+1,N+1);
for t=1:N
    phi=[u(t:-1:max(t-N,1));zeros(N+1-t,1)];
    P=(P-P*(phi*phi')*P/(lambda+phi'*P*phi))/lambda;
    K=P*phi;
    th(:,t+1)=th(:,t)+K*(y(t)-phi'*th(:,t));
    yhat(t)=phi'*th(:,t+1);
end
th=th(:,2:end)';
