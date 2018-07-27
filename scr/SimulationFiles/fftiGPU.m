function [dat] = fftiGPU(dat,NoConj,Conj)

dat=arrayfun(@leftToconj,dat,NoConj,Conj);
dat=ifft2(dat);
end

