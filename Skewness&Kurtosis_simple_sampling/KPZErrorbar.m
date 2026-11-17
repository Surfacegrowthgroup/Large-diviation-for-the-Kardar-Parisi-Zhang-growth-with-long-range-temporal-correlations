function [HT]=KPZErrorbar(N,thc)
%This function simulates the height data of the KPZ equation with long-range temporal correlations to calculate the skewness, kurtosis and corresponding error bar, using simple sampling.
%The input is the sample times (N), and the correlation exponent (thc). And the output is sample result (HT).
%The nonlinear coefficient and the simulation time T are included in the function.
L=512;
dt=0.01;
if thc==0
    T=9;
    nabla=19.5;
elseif thc==0.05
    T=27;
    nabla=14;
elseif thc==0.10
    T=30;
    nabla=15;
elseif thc==0.15
    T=27;
    nabla=9;
elseif thc==0.20
    T=28;
    nabla=13;
elseif thc==0.25
    T=30;
    nabla=8;
elseif thc==0.30
    T=35;
    nabla=10;
elseif thc==0.35
    T=26;
    nabla=7;
elseif thc==0.40
    T=95;
    nabla=10;
end
HT=zeros(N,L);
tic;
if thc==0
    parfor i=1:N
        theta1=normrnd(0,1,L,1,T*(1/dt));
        [A0]=KPZHblackbox(L,nabla,T,dt,theta1);
        HT(i,:)=A0;
    end
else
    parfor i=1:N
        theta1=rand(L,1,50,T*(1/dt));
        [A0]=KPZHblackbox(L,nabla,T,dt,(FFGN_revise(thc,T*(1/dt),theta1,L,1))/5);
        HT(i,:)=A0;
    end
end
toc;