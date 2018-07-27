function [dat] = ffti(dat,N1,N2,M1,M2)
% Reshape the data and make the inverse fast fourier transform of the data
dat=[dat; [conj(flipud(dat(2:N1,1))) conj(rot90(dat(2:N1,2:M2),2))]];
   
dat=ifft2(dat);

end

