function [dUdt] = GPU_dUdt_FullModel(yU,h0,M1,M2,N1,N2,kx,kz,delta,eta,zeta,k2cut,N,Conj,NoConj,GPU)
%function calculating the equation for the full isothermal model that will be solved by ode45 with the
%gpu calaulation . 
%
%For the meaning of the variables used please refer to the documentation  
%
%see also fft, ode45 
% 

yU(1)=h0+imag(yU(1))*1i;

%/* FOURIER SPACE: All the derivatives*/
      reconstru=reshape(yU,(N1+1)*M2,7);
datah = reconstru(:,1);
dataq = reconstru(:,2);
datap = reconstru(:,3);
dataS1 = reconstru(:,4);
dataS2 = reconstru(:,5);
dataR1 = reconstru(:,6);
dataR2 = reconstru(:,7);


clear reconstru

datah=reshape(datah,M2,N1+1)';
dataq=reshape(dataq,M2,N1+1)';
datap=reshape(datap,M2,N1+1)';
dataS1=reshape(dataS1,M2,N1+1)';
dataS2=reshape(dataS2,M2,N1+1)';
dataR1=reshape(dataR1,M2,N1+1)';
dataR2=reshape(dataR2,M2,N1+1)';

datah(1,1)=h0+imag(datah(1,1));


kx2=(kx.*kx);
kz2=(kz.*kz);
kxz=(kx.*kz);
k2=kx2+kz2;

M1M2=M1*M2;

if(GPU == true)
datah=gpuArray([datah; [(flipud(datah(2:N1,1))) (rot90(datah(2:N1,2:M2),2))]]);
dataq=gpuArray([dataq; [(flipud(dataq(2:N1,1))) (rot90(dataq(2:N1,2:M2),2))]]);
datap=gpuArray([datap; [(flipud(datap(2:N1,1))) (rot90(datap(2:N1,2:M2),2))]]);
dataS1=gpuArray([dataS1; [(flipud(dataS1(2:N1,1))) (rot90(dataS1(2:N1,2:M2),2))]]);
dataS2=gpuArray([dataS2; [(flipud(dataS2(2:N1,1))) (rot90(dataS2(2:N1,2:M2),2))]]);
dataR1=gpuArray([dataR1; [(flipud(dataR1(2:N1,1))) (rot90(dataR1(2:N1,2:M2),2))]]);
dataR2=gpuArray([dataR2; [(flipud(dataR2(2:N1,1))) (rot90(dataR2(2:N1,2:M2),2))]]);
else
datah=([datah; [(flipud(datah(2:N1,1))) (rot90(datah(2:N1,2:M2),2))]]);
dataq=([dataq; [(flipud(dataq(2:N1,1))) (rot90(dataq(2:N1,2:M2),2))]]);
datap=([datap; [(flipud(datap(2:N1,1))) (rot90(datap(2:N1,2:M2),2))]]);
dataS1=([dataS1; [(flipud(dataS1(2:N1,1))) (rot90(dataS1(2:N1,2:M2),2))]]);
dataS2=([dataS2; [(flipud(dataS2(2:N1,1))) (rot90(dataS2(2:N1,2:M2),2))]]);
dataR1=([dataR1; [(flipud(dataR1(2:N1,1))) (rot90(dataR1(2:N1,2:M2),2))]]);
dataR2=([dataR2; [(flipud(dataR2(2:N1,1))) (rot90(dataR2(2:N1,2:M2),2))]]);
end

      datahx=kx.*imag(datah)-kx.*real(datah)*1i;      
      dataqx=kx.*imag(dataq)-kx.*real(dataq)*1i;
      datapx=kx.*imag(datap)-kx.*real(datap)*1i;
      dataS1x=kx.*imag(dataS1)-kx.*real(dataS1)*1i;
      dataS2x=kx.*imag(dataS2)-kx.*real(dataS2)*1i;
      dataR1x=kx.*imag(dataR1)-kx.*real(dataR1)*1i;
      dataR2x=kx.*imag(dataR2)-kx.*real(dataR2)*1i;
      
      
      hx  = real(fftiGPU(datahx*M1M2,NoConj,Conj));
      px  = real(fftiGPU(datapx*M1M2,NoConj,Conj));
      qx  = real(fftiGPU(dataqx*M1M2,NoConj,Conj));
      s1x  = real(fftiGPU(dataS1x*M1M2,NoConj,Conj));
      s2x  = real(fftiGPU(dataS2x*M1M2,NoConj,Conj));
      r1x  = real(fftiGPU(dataR1x*M1M2,NoConj,Conj));
      r2x  = real(fftiGPU(dataR2x*M1M2,NoConj,Conj));
     
      %/* z-derivatives */
      datahz=kz.*imag(datah)-kz.*real(datah)*1i;
      dataqz=kz.*imag(dataq)-kz.*real(dataq)*1i;
      datapz=kz.*imag(datap)-kz.*real(datap)*1i;
      dataS1z=kz.*imag(dataS1)-kz.*real(dataS1)*1i;
      dataS2z=kz.*imag(dataS2)-kz.*real(dataS2)*1i;
      dataR1z=kz.*imag(dataR1)-kz.*real(dataR1)*1i;
      dataR2z=kz.*imag(dataR2)-kz.*real(dataR2)*1i;
      
      
      hZ  = real(fftiGPU(datahz*M1M2,NoConj,Conj));
      pz  = real(fftiGPU(datapz*M1M2,NoConj,Conj));
      qz  = real(fftiGPU(dataqz*M1M2,NoConj,Conj));
      s1z  = real(fftiGPU(dataS1z*M1M2,NoConj,Conj));
      s2z  = real(fftiGPU(dataS2z*M1M2,NoConj,Conj));
      r1z  = real(fftiGPU(dataR1z*M1M2,NoConj,Conj));
      r2z  = real(fftiGPU(dataR2z*M1M2,NoConj,Conj));
      
      %/* Double x-derivatives */
      datahxx=-kx2.*real(datah)-kx2.*imag(datah)*1i;
      dataqxx=-kx2.*real(dataq)-kx2.*imag(dataq)*1i;
      datapxx=-kx2.*real(datap)-kx2.*imag(datap)*1i;
      %dataS1xx=-kx2.*imag(dataS1)-kx2.*real(dataS1)*1i;
      %dataS2xx=-kx2.*imag(dataS2)-kx2.*real(dataS2)*1i;
      %dataR1xx=-kx2.*imag(dataR1)-kx2.*real(dataR1)*1i;
      %dataR2xx=-kx2.*imag(dataR2)-kx2.*real(dataR2)*1i;     
      
      hxx = real(fftiGPU(datahxx*M1M2,NoConj,Conj));
      pxx = real(fftiGPU(datapxx*M1M2,NoConj,Conj));
      qxx = real(fftiGPU(dataqxx*M1M2,NoConj,Conj));
      %s1xx  = real(fftiGPU(dataS1xx*M1M2,N1,N2,M1,M2,NoConj,Conj));
      %s2xx  = real(fftiGPU(dataS2xx*M1M2,N1,N2,M1,M2,NoConj,Conj));
      %r1xx  = real(fftiGPU(dataR1xx*M1M2,N1,N2,M1,M2,NoConj,Conj));
      %r2xx  = real(fftiGPU(dataR2xx*M1M2,N1,N2,M1,M2,NoConj,Conj));
      
      %/* Double z-derivatives */
      datahzz=-kz2.*real(datah)-kz2.*imag(datah)*1i;
      dataqzz=-kz2.*real(dataq)-kz2.*imag(dataq)*1i;
      datapzz=-kz2.*real(datap)-kz2.*imag(datap)*1i; 
      %dataS1zz=-kz2.*imag(dataS1)-kz2.*real(dataS1)*1i;
      %dataS2zz=-kz2.*imag(dataS2)-kz2.*real(dataS2)*1i;
      %dataR1zz=-kz2.*imag(dataR1)-kz2.*real(dataR1)*1i;
      %dataR2zz=-kz2.*imag(dataR2)-kz2.*real(dataR2)*1i;      
         
      hzz = real(fftiGPU(datahzz*M1M2,NoConj,Conj));
      pzz = real(fftiGPU(datapzz*M1M2,NoConj,Conj));
      qzz = real(fftiGPU(dataqzz*M1M2,NoConj,Conj));
      %s1zz  = real(fftiGPU(dataS1zz*M1M2,N1,N2,M1,M2,NoConj,Conj));
      %s2zz  = real(fftiGPU(dataS2zz*M1M2,N1,N2,M1,M2,NoConj,Conj));
      %r1zz  = real(fftiGPU(dataR1zz*M1M2,N1,N2,M1,M2,NoConj,Conj));
      %r2zz  = real(fftiGPU(dataR2zz*M1M2,N1,N2,M1,M2,NoConj,Conj));
      
      %/* Cross derivatives */  %CHECK VORZEICHEN
      datahxz=-kxz.*real(datah)-kxz.*imag(datah)*1i;
      dataqxz=-kxz.*real(dataq)-kxz.*imag(dataq)*1i;
      datapxz=-kxz.*real(datap)-kxz.*imag(datap)*1i;
      %dataS1xz=-kxz.*imag(dataS1)-kxz.*real(dataS1)*1i;
      %dataS2xz=-kxz.*imag(dataS2)-kxz.*real(dataS2)*1i;
      %dataR1xz=-kxz.*imag(dataR1)-kxz.*real(dataR1)*1i;
      %dataR2xz=-kxz.*imag(dataR2)-kxz.*real(dataR2)*1i;   
      
      hxz = real(fftiGPU(datahxz*M1M2,NoConj,Conj));
      pxz = real(fftiGPU(datapxz*M1M2,NoConj,Conj));
      qxz = real(fftiGPU(dataqxz*M1M2,NoConj,Conj));
      %s1xz  = real(fftiGPU(dataS1xz*M1M2,N1,N2,M1,M2,NoConj,Conj));
      %s2xz  = real(fftiGPU(dataS2xz*M1M2,N1,N2,M1,M2,NoConj,Conj));
      %r1xz  = real(fftiGPU(dataR1xz*M1M2,N1,N2,M1,M2,NoConj,Conj));
      %r2xz  = real(fftiGPU(dataR2xz*M1M2,N1,N2,M1,M2,NoConj,Conj));
      
      %/* Gradient of the laplacien of h (=curvature K) 
      datahLx=-(kx.*k2.*imag(datah)-kx.*k2.*real(datah)*1i);
      datahLz=-(kz.*k2.*imag(datah)-kz.*k2.*real(datah)*1i); 
      
       dhdt=arrayfun(@fU_dhdt,dataqx,datapz);
       
       
h   = real(fftiGPU(datah*M1M2,NoConj,Conj));
P   = real(fftiGPU(datap*M1M2,NoConj,Conj));
q   = real(fftiGPU(dataq*M1M2,NoConj,Conj));
s1  = real(fftiGPU(dataS1*M1M2,NoConj,Conj));
s2  = real(fftiGPU(dataS2*M1M2,NoConj,Conj));
r1  = real(fftiGPU(dataR1*M1M2,NoConj,Conj));
r2  = real(fftiGPU(dataR2*M1M2,NoConj,Conj));

datahLx = fftiGPU(datahLx*M1M2,NoConj,Conj);
datahLz = fftiGPU(datahLz*M1M2,NoConj,Conj);

% %/* REAL SPACE: */
h2=h.*h;
h3=h2.*h;

dqdt= (4.292582417582418e-6.*(520.*eta.*h2.*(1579.*pxz + 2027.*qxx + 448.*qzz) + ...
                6.*(65.*eta.*(2539.*hx.*hZ.*P + 3350.*hx.*hx.*q + 811.*hZ.*hZ.*q) - ...
                4160.*(27.*q + 308.*s1 + 1023.*s2) + ...
                1792.*delta.*(2.*hx.*q.*(13.*q - 26.*s1 - 21.*s2) + ...
                hZ.*(26.*P.*q - 26.*q.*r1 - 21.*q.*r2 - 26.*P.*s1 - 21.*P.*s2))) - ...
                3.*h.*(65.*eta.*(5447.*hxz.*P + 3515.*hZ.*px + 1583.*hx.*pz + 7234.*hxx.*q + ...
                1787.*hzz.*q + 6740.*hx.*qx + 1642.*hZ.*qz) + ...
                512.*delta.*(14.*pz.*(13.*q - 13.*s1 - 18.*s2) + ...
                q.*(364.*qx - 182.*r1z - 192.*r2z - 364.*s1x - 339.*s2x) + ...
                7.*(26.*P.*qz - 26.*qz.*r1 - 21.*qz.*r2 - 52.*qx.*s1 - 26.*P.*s1z - 57.*qx.*s2 - ...
                21.*P.*s2z))) - 224640.*h3.*(-1. + hx.*zeta) + ... 
                224640.*h3.*real(datahLx)))./(delta.*h2);

ds1dt= (0.000012487512487512488.*(-18018.*eta.*h2.*(pxz + qxx) + ...
                3.*(1001.*eta.*(35.*hx.*hZ.*P + 62.*hx.*hx.*q + 27.*hZ.*hZ.*q) + ...
                16.*delta.*(hZ.*(9.*q.*(182.*r1 - 93.*r2) + P.*(-143.*q + 1638.*s1 - 837.*s2)) + ...
                hx.*q.*(-143.*q + 3276.*s1 - 1674.*s2)) - 8008.*(q + 84.*(s1 + s2))) + ...
                h.*(3003.*eta.*(9.*hxz.*P - 6.*hZ.*px - 21.*hx.*pz + 7.*hxx.*q - 2.*hzz.*q - ...
                46.*hx.*qx - 19.*hZ.*qz) + 8.*delta.*(-2.*pz.*(286.*q + 4459.*s1 - 7146.*s2) + ...
                q.*(286.*qx + 2730.*r1z + 5535.*r2z - 7098.*s1x + 10557.*s2x) + ...
                2.*(429.*P.*qz - 4914.*qz.*r1 + 2511.*qz.*r2 - 9373.*qx.*s1 - 4914.*P.*s1z + ...
                9657.*qx.*s2 + 2511.*P.*s2z))) - 8008.*h3.*(-1. + hx.*zeta) +   ...
                8008.*h3.*real(datahLx)))./(delta.*h2) ;

ds2dt= (-1.691017316017316e-6.*(-147576.*eta.*h2.*(pxz + qxx) + ...
                6.*(143.*eta.*(215.*hx.*hZ.*P + 494.*hx.*hx.*q + 279.*hZ.*hZ.*q) + ...
                9152.*(q + 84.*s1 + 909.*s2) + 8960.*delta.*(2.*hx.*q.*(2.*s1 - 9.*s2) + ...
                hZ.*(2.*q.*r1 - 9.*q.*r2 + 2.*P.*s1 - 9.*P.*s2))) + ...
                h.*(429.*eta.*(461.*hxz.*P + 41.*hZ.*px - 379.*hx.*pz + 438.*hxx.*q - 23.*hzz.*q - ...
                804.*hx.*qx - 466.*hZ.*qz) + 512.*delta.*(105.*qz.*(-2.*r1 + 9.*r2) + ...
                3.*q.*(28.*r1z - 27.*r2z - 42.*s1x + 288.*s2x) + ...
                35.*(8.*pz.*s1 + 2.*qx.*s1 - 6.*P.*s1z + 30.*pz.*s2 + 57.*qx.*s2 + 27.*P.*s2z))) + ...
                18304.*h3.*(-1. + hx.*zeta) - 18304.*h3.*real(datahLx)))./(delta.*h2) ;
  
            
dpdt= (4.292582417582418e-6.*(279552.*delta.*hZ.*P.*P + ...
                65.*eta.*(15234.*hx.*hZ.*q + 8.*h2.*(448.*pxx + 2027.*pzz + 1579.*qxz) - ...
                3.*h.*(1642.*hx.*px + 6740.*hZ.*pz + 5447.*hxz.*q + 1583.*hZ.*qx + 3515.*hx.*qz)) + ...
                3.*P.*(65.*eta*(1622.*hx.*hx - 1787.*h.*hxx + 6700.*hZ.*hZ - 7234.*h.*hzz) - 128.*(1755. + ...
                4.*delta.*(7.*(-26.*hx.*q + 52.*hZ.*r1 + 42.*hZ.*r2 + 26.*hx.*s1 + 21.*hx.*s2) + ...
                h.*(364.*pz + 182.*qx - 364.*r1z - 339.*r2z - 182.*s1x - 192.*s2x)))) - ...
                384.*(28.*delta.*(hx.*q.*(26.*r1 + 21.*r2) + ...
                h.*(26.*px.*q - 52.*pz.*r1 - 26.*qx.*r1 - 26.*q.*r1x - 57.*pz.*r2 - 36.*qx.*r2 - ...
                21.*q.*r2x - 26.*px.*s1 - 21.*px.*s2)) + 65.*(308.*r1 + 1023.*r2 + 9.*h3.*hZ.*zeta )) + ...
                224640.*h3.*real(datahLz)))./(delta*h2);        
   
dr1dt= (0.000012487512487512488*(-6864.*delta.*hZ.*P.*P - ...
                3003.*eta.*(-35.*hx.*hZ.*q + 6.*h2.*(pzz + qxz) + ...
                h.*(19.*hx.*px + 46.*hZ.*pz - 9.*hxz.*q + 21.*hZ.*qx + 6.*hx.*qz)) - ...
                2.018016e6.*(r1 + r2) + 16.*delta.*(27.*hx.*q.*(182.*r1 - 93.*r2) + ...
                h.*(-9373.*pz.*r1 - 4459.*qx.*r1 - 4914.*q.*r1x + 9657.*pz.*r2 + 7146.*qx.*r2 + ...
                2511.*q.*r2x + 3.*px.*(143.*q - 1638.*s1 + 837.*s2))) + ...
                P.*(-24024. + 3003.*eta.*(27.*hx.*hx - 2.*h.*hxx + 62.*hZ.*hZ + ...
                7.*h.*hzz) + 8.*delta.*(6.* ...
                (18.*hZ.*(182.*r1 - 93.*r2) + hx.*(-143.*q + 1638.*s1 - 837.*s2)) + ...
                h.*(286.*pz - 572.*qx - 7098.*r1z + 10557.*r2z + 2730.*s1x + 5535.*s2x))) - ...
                8008.*h3.*hZ.*zeta + 8008.*h3.*real(datahLz)))./(delta.*h2) ;
            
dr2dt= (-1.691017316017316e-6.*(-429.*eta.*(-430.*hx.*hZ.*q + 344.*h2.*(pzz + qxz) + ...
                h.*(466.*hx.*px + 804.*hZ.*pz - 461.*hxz.*q + 379.*hZ.*qx - 41.*hx.*qz)) + ...
                164736.*(28.*r1 + 303.*r2) + 17920.*delta.*(3.*hx.*q.*(2.*r1 - 9.*r2) + ...
                h.*(2.*pz.*r1 + 8.*qx.*r1 - 6.*q.*r1x + 57.*pz.*r2 + 30.*qx.*r2 + 27.*q.*r2x - ...
                6.*px.*s1 + 27.*px.*s2)) + ...
                3.*P.*(18304. + 143.*eta.*(558.*hx.*hx - 23.*h.*hxx + 988.*hZ.*hZ + ...
                438.*h.*hzz) + 512.*delta.*(70.*hZ.*(2.*r1 - 9.*r2) + 35.*hx.*(2.*s1 - 9.*s2) + ...
                h.*(-42.*r1z + 288.*r2z + 28.*s1x - 27.*s2x))) + ...
                18304.*h3.*hZ.*zeta - 18304.*h3.*real(datahLz)))./(delta.*h2) ;
              
% /* FOURIER SPACE (attention /(M1*M2)!  */
dqdt=fft2(dqdt);
dpdt=fft2(dpdt);
ds1dt=fft2(ds1dt);
ds2dt=fft2(ds2dt);
dr1dt=fft2(dr1dt);
dr2dt=fft2(dr2dt);

%    /*cut to avoid aliasing errors*/
dqdt=dqdt.*k2cut/M1M2;
dpdt=dpdt.*k2cut/M1M2;
ds1dt=ds1dt.*k2cut/M1M2;
ds2dt=ds2dt.*k2cut/M1M2;
dr1dt=dr1dt.*k2cut/M1M2;
dr2dt=dr2dt.*k2cut/M1M2;

dhdt=conj(dhdt);
dqdt=conj(dqdt);
dpdt=conj(dpdt);
ds1dt=conj(ds1dt);
ds2dt=conj(ds2dt);
dr1dt=conj(dr1dt);
dr2dt=conj(dr2dt);


dpdt=gather(dpdt);
dhdt=gather(dhdt); 
dqdt=gather(dqdt);
ds1dt=gather(ds1dt);
ds2dt=gather(ds2dt);
dr1dt=gather(dr1dt);
dr2dt=gather(dr2dt);
  

dhdt=reshape(dhdt',1,[]);
dqdt=reshape(dqdt',1,[]);
dpdt=reshape(dpdt',1,[]);
ds1dt=reshape(ds1dt',1,[]);
ds2dt=reshape(ds2dt',1,[]);
dr1dt=reshape(dr1dt',1,[]);
dr2dt=reshape(dr2dt',1,[]);

dUdt=[dhdt(1:N/2) dqdt(1:N/2) dpdt(1:N/2) ds1dt(1:N/2) ds2dt(1:N/2) dr1dt(1:N/2) dr2dt(1:N/2)]';
end