function [HT,acceptence,accwre]=KPZHcal(thd,N,thc)
%This function is for simulating the height data of the long-range temporal correlation KPZ equation using the large deviation algorithm.
%The input is the biased parameter (thd), the sample times (N), and the correlated exponent (thc). And the output is sample result (HT), the acceptance rate in the sampling process (acceptence), and the acceptance rate in the test sampling process (accwre).
L=512;
if thc==0
    T=9;
    nabla=19.5;
elseif thc==0.05
    T=29;
    nabla=14;
elseif thc==0.15
    T=27;
    nabla=9;
elseif thc==0.25
    T=30;
    nabla=8;
elseif thc==0.35
    T=26;
    nabla=7;
elseif thc==0.40
    T=95;
    nabla=7;
end
dt=0.01;
acc=0;
if thc~=0
    if thd~=0
        [P,accwre]=dichotomyKPZH(0.02,thd,thc,nabla);
        HT=zeros(N,1);
        n=round(P*L*T*(1/dt)*50);
        tic;
        theta0=rand(L,1,50,T*(1/dt));
        [A0]=KPZHblackbox(L,nabla,T,dt,(FFGN_revise(thc,T*(1/dt),theta0,L,1))/5);
        for nn=1:N
            theta1=theta0;
            theta1(randperm(L*T*(1/dt)*50,n))=rand(n,1);
            [AT]=KPZHblackbox(L,nabla,T,dt,(FFGN_revise(thc,T*(1/dt),theta1,L,1))/5);
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
        accwre=1;
        HT=zeros(N,1);
        tic;
        theta0=rand(L,1,50,T*(1/dt));
        [A0]=KPZHblackbox(L,nabla,T,dt,(FFGN_revise(thc,T*(1/dt),theta0,L,1))/5);
        for nn=1:N
            theta1=rand(L,1,50,T*(1/dt));
            [AT]=KPZHblackbox(L,nabla,T,dt,(FFGN_revise(thc,T*(1/dt),theta1,L,1))/5);
            hW=min(1,exp(-thd*(AT-A0)));
            if rand(1,1)<=hW
                % theta0=theta1;
                A0=AT;
                acc=acc+1;
            end
            HT(nn)=A0;
        end
        toc;
    end
elseif thc==0
    if thd~=0
        [P,accwre]=dichotomyKPZH(0.01,thd,thc,nabla);
        HT=zeros(N,1);
        n=round(P*L*T*(1/dt));
        tic;
        theta0=normrnd(0,1,L,1,T*(1/dt));
        [A0]=KPZHblackbox(L,nabla,T,dt,theta0);
        for nn=1:N
            theta1=theta0;
            theta1(randperm(L*T*(1/dt),n))=normrnd(0,1,n,1);
            [AT]=KPZHblackbox(L,nabla,T,dt,theta1);
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
        accwre=1;
        HT=zeros(N,1);
        tic;
        theta0=normrnd(0,1,L,1,T*(1/dt));
        [A0]=KPZHblackbox(L,nabla,T,dt,theta0);
        for nn=1:N
            theta1=normrnd(0,1,L,1,T*(1/dt));
            [AT]=KPZHblackbox(L,nabla,T,dt,theta1);
            hW=min(1,exp(-thd*(AT-A0)));
            if rand(1,1)<=hW
                % theta0=theta1;
                A0=AT;
                acc=acc+1;
            end
            HT(nn)=A0;
        end
        toc;
    end
end
acceptence=acc/N;
end