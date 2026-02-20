function [noise_result]=FFGN_revise(thc,T,thrn,sz1,sz2)
%This function is for generating the long-range temporally correlated noise using the fast fractional Gaussian noise (FFGN) method.
%The input is the correlated exponent (th), the length of correlated sequence (N), the all-needed random number (thrn), the system size in d-dimensional substrate (sz1, sz2), and the output is the long-range temporal correlated noise with the selected substrate (noise_result).
%The matrix size is thrn (sz1*sz2*ga*t), noise_result (sz1*sz2*T).
%ga=50;
noise_result=zeros(sz1,sz2,T);
un=6*2.^(-(1:50));
rn=exp(-1*un);
Wn=(12*(1-rn.^2)*(2^(0.5-thc)-2^(thc-0.5)).*(un.^(1-2*thc)))/(gamma(2-2*thc));
Wn=permute(sqrt(Wn),[1,3,2]);
rn=permute(rn,[1,3,2]);
X0=((1-(rn.^2)).^(-0.5)).*((thrn(:,:,:,1)-0.5));
noise_result(:,:,1)=sum(Wn.*X0,3);
for i=2:T
    X1=rn.*X0+((thrn(:,:,:,i)-0.5));
    noise_result(:,:,i)=sum(Wn.*X1,3);
    X0=X1;
end
end

