function [acc]=acceptenceKPZH(P,thd,thc,nabla)
%The function is used to test the acceptance rate during the Markov chain with the input parameter condition. 
%The replacement rate of random numbers (P), the biased parameter (thd), the correlated exponent (thc), and the nonlinear parameter (nabla) is inputing. And the function can output the acceptance rate (acc).
dt=0.01;
if thc==0
    T=9;
elseif thc==0.05
    T=29;
elseif thc==0.15
    T=27;
elseif thc==0.25
    T=30;
elseif thc==0.35
    T=26;
elseif thc==0.38
    T=37;
elseif thc==0.40
    T=95;
end
N=10000;
L=512;
if thc~=0
    n=round(P*L*T*(1/dt)*50);
    theta0=rand(L,1,50,T*(1/dt));
    [A0]=KPZHblackbox(L,nabla,T,dt,(FFGN_revise(thc,T*(1/dt),theta0,L,1))/5);
    acce=0;
    raa=rand(N,1);
    for nn=1:N
        theta1=theta0;
        theta1(randperm(L*T*(1/dt)*50,n))=rand(n,1);
        [AT]=KPZHblackbox(L,nabla,T,dt,(FFGN_revise(thc,T*(1/dt),theta1,L,1))/5);
        hW=min(1,exp(-thd*(AT-A0)));
        if  raa(nn)<=hW
            theta0=theta1;
            A0=AT;
            acce=acce+1;
        end
    end
else
    n=round(P*L*T*(1/dt));
    theta0=normrnd(0,1,L,1,T*(1/dt));
    [A0]=KPZHblackbox(L,nabla,T,dt,theta0);
    acce=0;
    raa=rand(N,1);
    for nn=1:N
        theta1=theta0;
        theta1(randperm(L*T*(1/dt),n))=normrnd(0,1,n,1);
        [AT]=KPZHblackbox(L,nabla,T,dt,theta1);
        hW=min(1,exp(-thd*(AT-A0)));
        if  raa(nn)<=hW
            theta0=theta1;
            A0=AT;
            acce=acce+1;
        end
    end
end

acc=acce/N;
end