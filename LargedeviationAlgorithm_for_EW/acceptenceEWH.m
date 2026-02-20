function [acc]=acceptenceEWH(P,thd,thc)
%The function is used to test the acceptance rate during the Markov chain with the input parameter condition. 
%The replacement rate of random numbers (P), the biased parameter (thd), the correlated exponent (thc) is inputing. And the function can output the acceptance rate (acc).
dt=0.05;
T=10;
N=10000;
L=512;
n=round(P*512*T*(1/dt)*50);
theta0=rand(L,1,50,T*(1/dt));
[A0]=EWHblackbox(L,T,dt,FFGN_revise(thc,T*(1/dt),theta0,L,1));
acce=0;
raa=rand(N,1);
for nn=1:N
    theta1=theta0;
    theta1(randperm(L*T*(1/dt)*50,n))=rand(n,1);
    [AT]=EWHblackbox(L,T,dt,FFGN_revise(thc,T*(1/dt),theta1,L,1));
    hW=min(1,exp(-thd*(AT-A0)));
    if  raa(nn)<=hW
        theta0=theta1;
        A0=AT;
        acce=acce+1;
    end
end
acc=acce/N;
end