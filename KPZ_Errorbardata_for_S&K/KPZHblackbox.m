function [A]=KPZHblackbox(L,nabla,T,dt,noise)
%The function is the calculation process of the KPZ equation in the presence of long-range temporal correlations. 
%The input is the system size (L), the nonlinear parameter (nabla), the total time step (T), the time interval (dt), the noise data of all iteration steps (noise), and the output is the height result of all position in time T (H). 
A=zeros(L,1);
for i=1:T
    for k=1:(1/dt)
        [A]=KPZiteration(nabla,A,noise(:,:,(i-1)*(1/dt)+k),dt);
    end
end
end