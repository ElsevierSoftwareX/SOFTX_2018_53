function [datah,dataq,datap]   = reconstruct_U(U,N1,N2,M1,M2)
%Extracting Data from U Vector
sizeU=size(U(end,:));

if sizeU(2) == M2*(N1+1)*7
    reconstr=reshape(U(end,:),(N1+1)*M2,7);  %Full second-order model
else
    reconstr=reshape(U(end,:),(N1+1)*M2,3);  %Simplified model
end

datah=reshape(reconstr(:,1),M2,N1+1)';
dataq=reshape(reconstr(:,2),M2,N1+1)';
datap=reshape(reconstr(:,3),M2,N1+1)';

%Reconstruct h,q,p by inverse fourier transformation
datah   = ffti(datah*M1*M2,N1,N2,M1,M2);
dataq   = ffti(dataq*M1*M2,N1,N2,M1,M2);
datap   = ffti(datap*M1*M2,N1,N2,M1,M2);
end

