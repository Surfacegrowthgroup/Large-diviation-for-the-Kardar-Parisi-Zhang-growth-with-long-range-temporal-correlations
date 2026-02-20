function [bestP,acc11]=dichotomyEWH(dA,thd,thc)
%The function is for finding the best replacement rate of random numbers to ensure the acceptance rate is close to 0.5. This algorithm is based on the bisection method.
%The input is the acceptable deviation (dA), biased parameter (thd), correlated exponent (thc). The best replacement rate (bestP) and the corresponding acceptance rate (acc11) is outputing in the function.
acc11=0;
P1=0;
P=0.5;
P0=P;
while abs((acc11-0.5))>=dA
    acc11=acceptenceEWH(P,thd,thc);
    bestP=P;
if (acc11-0.5)<0
    P=P0-(abs(P0-P1))/2;
    P1=P0;
    P0=P;
elseif (acc11-0.5)>0
    P=P0+(abs(P0-P1))/2;
    P1=P0;
    P0=P;
end

end
end