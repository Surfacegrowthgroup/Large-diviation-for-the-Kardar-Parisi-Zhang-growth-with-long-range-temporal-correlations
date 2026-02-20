function [At]=EWiteration(A0,noise,dt)
%The function is for the iteration of the EW equation driven by the temporally correlated noise. 
%The input is the original height in time t (A0), noise data (noise), time iteration (dt), and the output is the height in the next time t+1 (At). 
L=length(A0);
At=((A0([2:L,1])-2*A0(1:L)+A0([L,1:(L-1)])))*dt+A0+noise*dt;
end