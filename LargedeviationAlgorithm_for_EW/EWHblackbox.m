function [A0]=EWHblackbox(L,T,dt,noise)
%The function is the calculation process of the EW equation in the presence of long-range temporal correlations. 
%The input is the system size (L), the total time step (T), the time interval (dt), the noise data of all iteration steps (noise), and the output is the height result of one point in time T (A0). 
A=zeros(L,1);
for i=1:T
    for k=1:(1/dt)
        [A]=EWiteration(A,noise(:,:,(i-1)*(1/dt)+k),dt);
    end
end
A0=A(256);
end

