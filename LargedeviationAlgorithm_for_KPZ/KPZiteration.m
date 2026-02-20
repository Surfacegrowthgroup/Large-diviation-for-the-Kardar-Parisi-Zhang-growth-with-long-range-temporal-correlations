function [A1]=KPZiteration(nabla,A0,noise,dt)
%The function is for the iteration of the KPZ equation driven by the temporally correlated noise. 
%The input is the nonlinear parameter (nabla), original height in time t (A0), noise data (noise), time iteration (dt), and the output is the height in the next time t+1 (A1). 
L=length(A0);
linear=(A0([2:L,1])-2*A0(1:L)+A0([L,1:L-1]));
nonlinear=((A0([2:L,1])-A0(1:L)).^2+(A0([2:L,1])-A0(1:L)).*(A0(1:L)-A0([L,1:L-1]))+(A0(1:L)-A0([L,1:L-1])).^2)/3;
A1=linear*dt+nabla*nonlinear*dt+A0(1:L)+noise(1:L)*dt;
end