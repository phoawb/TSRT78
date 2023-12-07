function [th, yhat] = rekid_lms(y,u,N,mu)

ny=length(y);
th=zeros(N+1,ny+1); yhat=zeros(ny,1);
for t=1:ny
    phi=[u(t:-1:max(t-N,1));zeros(N+1-t,1)];
    K=mu*phi;
    th(:,t+1)=th(:,t)+K*(y(t)-phi'*th(:,t));
    yhat(t)=phi'*th(:,t+1);
end
th=th(:,2:end)';
