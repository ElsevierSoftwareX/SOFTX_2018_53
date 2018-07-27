function dat = leftToconj(dat,NoConj,Conj)

dat=dat.*NoConj+conj(dat.*Conj);
end