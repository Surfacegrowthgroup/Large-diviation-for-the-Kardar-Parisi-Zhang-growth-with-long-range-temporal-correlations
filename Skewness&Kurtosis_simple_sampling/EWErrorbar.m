function [HT]=EWErrorbar(N,thc)
%This function simulates the height data of the EW equation with long-range temporal correlations to calculate the skewness, kurtosis and corresponding error bar, using simple sampling.
%The input is the sample times (N), and the correlation exponent (thc). And the output is sample result (HT).
L=512;
T=10;
dt=0.05;
HT=zeros(N,L);
tic;
if thc==0
    parfor i=1:N
        theta1=normrnd(0,1,L,1,T*(1/dt))
        [AT]=EWHblackbox(L,T,dt,theta1);
        HT(i,:)=AT;
    end
else
    parfor i=1:N
        theta1=rand(L,1,50,T*(1/dt));
        [AT]=EWHblackbox(L,T,dt,FFGN_revise(thc,T*(1/dt),theta1,L,1));
        HT(i,:)=AT;
    end
end
toc;