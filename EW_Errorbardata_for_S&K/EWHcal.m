function [HESK]=EWHcal(N,thc)
%This function is for simulating the height data of the long-range temporal correlation EW equation for estimate the error bars of skewness and kurtosis.
%The input is the repeat times (N), the correlated exponent (thc), and the output is the sample result of height (HESK).
L=512;
T=10;
dt=0.05;
HESK=zeros(N,L);
if thc==0
    parfor i=1:N
        theta1=normrnd(0,1,L,1,T*(1/dt))
        [AT]=EWHblackbox(L,T,dt,theta1);
        HESK(i,:)=AT;
    end
else
    parfor i=1:N
        theta1=rand(L,1,50,T*(1/dt));
        [AT]=EWHblackbox(L,T,dt,FFGN_revise(thc,T*(1/dt),theta1,L,1));
        HESK(i,:)=AT;
    end
end