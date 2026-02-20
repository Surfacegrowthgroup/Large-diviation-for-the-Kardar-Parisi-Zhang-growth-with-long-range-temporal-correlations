function [HT,acceptence,accwre]=EWHcal(thd,N,thc)
%This function is for simulating the height data of the long-range temporal correlation EW equation using the large deviation algorithm.
%The input is the biased parameter (thd), the sample times (N), and the correlated exponent (thc). And the output is sample result (HT), the acceptance rate in the sampling process (acceptence), and the acceptance rate in the test sampling process (accwre).
L=512;
T=10;
dt=0.05;
if thd~=0
    [P,accwre]=dichotomyEWH(0.02,thd,thc);
    HT=zeros(N,1);
    acc=0;
    n=round(P*L*T*(1/dt)*50);
    tic;
    theta0=rand(L,1,50,T*(1/dt));
    [A0]=EWHblackbox(L,T,dt,FFGN_revise(thc,T*(1/dt),theta0,L,1));
    for nn=1:N
        theta1=theta0;
        theta1(randperm(L*T*(1/dt)*50,n))=rand(n,1);
        [AT]=EWHblackbox(L,T,dt,FFGN_revise(thc,T*(1/dt),theta1,L,1));
        hW=min(1,exp(-thd*(AT-A0)));
        if rand(1,1)<=hW
            theta0=theta1;
            A0=AT;
            acc=acc+1;
        end
        HT(nn)=A0;
    end
    toc;


else
    P=1;
    accwre=0;
    HT=zeros(N,1);
    acc=0;
    tic;
    theta0=rand(L,1,50,T*(1/dt));
    [A0]=EWHblackbox(L,T,dt,FFGN_revise(thc,T*(1/dt),theta0,L,1));
    for nn=1:N
        theta1=rand(L,1,50,T*(1/dt));
        [AT]=EWHblackbox(L,T,dt,FFGN_revise(thc,T*(1/dt),theta1,L,1));
        hW=min(1,exp(-thd*(AT-A0)));
        if rand(1,1)<=hW
            theta0=theta1;
            A0=AT;
            acc=acc+1;
        end
        HT(nn)=A0;
    end
    toc;
end
acceptence=acc/N;