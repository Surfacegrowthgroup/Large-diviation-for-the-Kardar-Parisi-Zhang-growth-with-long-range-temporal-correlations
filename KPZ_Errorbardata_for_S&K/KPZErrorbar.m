function [HSK]=KPZErrorbar(N,thc)
%This function is for simulating the height data of the long-range temporal correlation KPZ equation for estimate the error bars of skewness and kurtosis.
%The input is the repeat times (N), the correlated exponent (thc), and the output is the sample result of height (HSK).
L=512;
dt=0.01;
if thc==0
    T=9;
    nabla=19.5;
elseif thc==0.05
    T=29;
    nabla=14;
elseif thc==0.10
    T=30;
    nabla=19.5;
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
    nabla=7;
end
HSK=zeros(N,L);
if thc==0
    parfor i=1:N
        theta1=normrnd(0,1,L,1,T*(1/dt));
        [A0]=KPZHblackbox(L,nabla,T,dt,theta1);
        HSK(i,:)=A0;
    end
else
    parfor i=1:N
        theta1=rand(L,1,50,T*(1/dt));
        [A0]=KPZHblackbox(L,nabla,T,dt,(FFGN_revise(thc,T*(1/dt),theta1,L,1))/5);
        HSK(i,:)=A0;
    end
end