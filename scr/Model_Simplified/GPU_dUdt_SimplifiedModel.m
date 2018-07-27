function [dUdt] = GPU_dUdt_SimplifiedModel(yU,h0,M1,M2,N1,N2,kx,kz,delta,eta,zeta,k2cut,N,Conj,NoConj,GPU)
%function calculating the equation for the simplified model that will be solved by ode45
% This variant will use GPU calculation 
%
yU(1)=h0+imag(yU(1))*1i;

%/* FOURIER SPACE: All the derivatives*/
      A=reshape(yU,(N1+1)*M2,3);
datah = A(:,1);
dataq = A(:,2);
datap = A(:,3);

datah=reshape(datah,M2,N1+1)';
dataq=reshape(dataq,M2,N1+1)';
datap=reshape(datap,M2,N1+1)';

kx2=(kx.*kx);
kz2=(kz.*kz);
kxz=(kx.*kz);
k2=kx2+kz2;

M1M2=M1*M2;

datah=gpuArray([datah; [(flipud(datah(2:N1,1))) (rot90(datah(2:N1,2:M2),2))]]);
dataq=gpuArray([dataq; [(flipud(dataq(2:N1,1))) (rot90(dataq(2:N1,2:M2),2))]]);
datap=gpuArray([datap; [(flipud(datap(2:N1,1))) (rot90(datap(2:N1,2:M2),2))]]);

      datahx=kx.*imag(datah)-kx.*real(datah)*1i;      
      dataqx=kx.*imag(dataq)-kx.*real(dataq)*1i;
      datapx=kx.*imag(datap)-kx.*real(datap)*1i;     
     
      hx  = real(fftiGPU(datahx*M1M2,NoConj,Conj));
      px  = real(fftiGPU(datapx*M1M2,NoConj,Conj));
      qx  = real(fftiGPU(dataqx*M1M2,NoConj,Conj));
      
      %/* z-derivatives */
      datahz=kz.*imag(datah)-kz.*real(datah)*1i;
      dataqz=kz.*imag(dataq)-kz.*real(dataq)*1i;
      datapz=kz.*imag(datap)-kz.*real(datap)*1i;
 
      hZ  = real(fftiGPU(datahz*M1M2,NoConj,Conj));
      pz  = real(fftiGPU(datapz*M1M2,NoConj,Conj));
      qz  = real(fftiGPU(dataqz*M1M2,NoConj,Conj));
      
      %/* Double x-derivatives */
      datahxx=-kx2.*real(datah)-kx2.*imag(datah)*1i;
      dataqxx=-kx2.*real(dataq)-kx2.*imag(dataq)*1i;
      datapxx=-kx2.*real(datap)-kx2.*imag(datap)*1i;
 
      hxx = real(fftiGPU(datahxx*M1M2,NoConj,Conj));
      pxx = real(fftiGPU(datapxx*M1M2,NoConj,Conj));
      qxx = real(fftiGPU(dataqxx*M1M2,NoConj,Conj));
      
      %/* Double z-derivatives */
      datahzz=-kz2.*real(datah)-kz2.*imag(datah)*1i;
      dataqzz=-kz2.*real(dataq)-kz2.*imag(dataq)*1i;
      datapzz=-kz2.*real(datap)-kz2.*imag(datap)*1i; 
      
      hzz = real(fftiGPU(datahzz*M1M2,NoConj,Conj));
      pzz = real(fftiGPU(datapzz*M1M2,NoConj,Conj));
      qzz = real(fftiGPU(dataqzz*M1M2,NoConj,Conj));
      
      %/* Cross derivatives */  %CHECK VORZEICHEN
      datahxz=-kxz.*real(datah)-kxz.*imag(datah)*1i;
      dataqxz=-kxz.*real(dataq)-kxz.*imag(dataq)*1i;
      datapxz=-kxz.*real(datap)-kxz.*imag(datap)*1i;
      
      hxz = real(fftiGPU(datahxz*M1M2,NoConj,Conj));
      pxz = real(fftiGPU(datapxz*M1M2,NoConj,Conj));
      qxz = real(fftiGPU(dataqxz*M1M2,NoConj,Conj));
 
      %/* Gradient of the laplacien of h (=curvature K) 
      datahLx=-(kx.*k2.*imag(datah)-kx.*k2.*real(datah)*1i);
      datahLz=-(kz.*k2.*imag(datah)-kz.*k2.*real(datah)*1i); 
      
       dhdt=arrayfun(@GPU_dhdt,dataqx,datapz);
       
      
h   = real(fftiGPU(datah*M1M2,NoConj,Conj));
P   = real(fftiGPU(datap*M1M2,NoConj,Conj));
q   = real(fftiGPU(dataq*M1M2,NoConj,Conj));

datahLx = fftiGPU(datahLx*M1M2,NoConj,Conj);
datahLz = fftiGPU(datahLz*M1M2,NoConj,Conj);

% %/* REAL SPACE: */

h2=h.*h;
h3=h2.*h;

corr=1.0;  %/*coeff to reflect corrections inertial order 2*/

dqdt=arrayfun(@GPU_dqdt,h,h2,hx,hZ,hxz,hxx,hzz,q,qx,qz,qxx,qzz,P,pz,px,pxz,delta,eta,zeta,datahLx,corr);
dpdt=arrayfun(@GPU_dpdt,h,h2,h3,hx,hZ,hxz,hxx,hzz,q,qx,qz,qxz,P,pz,px,pxx,pzz,datahLz,delta,eta,zeta);


% /* FOURIER SPACE (attention /(M1*M2)!  */
dqdt=fft2(dqdt);
dpdt=fft2(dpdt);

%    /* cut to avoid aliasing errors  */
dqdt=dqdt.*k2cut/M1M2;
dpdt=dpdt.*k2cut/M1M2;

dhdt=conj(dhdt);
dqdt=conj(dqdt);
dpdt=conj(dpdt);

if(GPU == true)
dpdt=gather(dpdt);
dhdt=gather(dhdt); 
dqdt=gather(dqdt);
end

dhdt=reshape(dhdt',1,[]);
dqdt=reshape(dqdt',1,[]);
dpdt=reshape(dpdt',1,[]);

dUdt=[dhdt(1:N/2) dqdt(1:N/2) dpdt(1:N/2)]';

end