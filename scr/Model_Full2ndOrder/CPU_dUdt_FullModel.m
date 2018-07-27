%This function evaluates the derivative of the vector U based on the
%full second-order model that will be solved by ode45
%for the meaning of the variable used please refer to the documentation  
%
%see also fft, ode45
% 

function [dUdt] = CPU_dUdt_FullModel(yU,h0,M1,M2,N1,N2,kx,kz,delta,eta,zeta,k2cut,N)


%/* FOURIER SPACE: All the derivatives*/
%First, all variables are reconstructed from the state vector U
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

datah(1,1)=h0+imag(datah(1,1)); %Fixes the average film thicknes to the given value

%Cut-off matrix for the aliasing filter
kx2=kx.*kx;
kz2=kz.*kz;
kxz=kx.*kz;
k2=kx2+kz2;

M1M2=M1*M2;

%The following lines calculate the derivatives in all spatial directions
%and convertes the solution in real space

      %/* x-derivatives */
      datahx=kx.*imag(datah)-kx.*real(datah)*1i;      
      dataqx=kx.*imag(dataq)-kx.*real(dataq)*1i;
      datapx=kx.*imag(datap)-kx.*real(datap)*1i;
      dataS1x=kx.*imag(dataS1)-kx.*real(dataS1)*1i;
      dataS2x=kx.*imag(dataS2)-kx.*real(dataS2)*1i;
      dataR1x=kx.*imag(dataR1)-kx.*real(dataR1)*1i;
      dataR2x=kx.*imag(dataR2)-kx.*real(dataR2)*1i;
      
      hx  = real(ffti(datahx*M1M2,N1,N2,M1,M2));
      px  = real(ffti(datapx*M1M2,N1,N2,M1,M2));
      qx  = real(ffti(dataqx*M1M2,N1,N2,M1,M2));
      s1x  = real(ffti(dataS1x*M1M2,N1,N2,M1,M2));
      s2x  = real(ffti(dataS2x*M1M2,N1,N2,M1,M2));
      r1x  = real(ffti(dataR1x*M1M2,N1,N2,M1,M2));
      r2x  = real(ffti(dataR2x*M1M2,N1,N2,M1,M2));
      
      
      %/* z-derivatives */
      datahz=kz.*imag(datah)-kz.*real(datah)*1i;
      dataqz=kz.*imag(dataq)-kz.*real(dataq)*1i;
      datapz=kz.*imag(datap)-kz.*real(datap)*1i;
      dataS1z=kz.*imag(dataS1)-kz.*real(dataS1)*1i;
      dataS2z=kz.*imag(dataS2)-kz.*real(dataS2)*1i;
      dataR1z=kz.*imag(dataR1)-kz.*real(dataR1)*1i;
      dataR2z=kz.*imag(dataR2)-kz.*real(dataR2)*1i;
      
      hZ  = real(ffti(datahz*M1M2,N1,N2,M1,M2));
      pz  = real(ffti(datapz*M1M2,N1,N2,M1,M2));
      qz  = real(ffti(dataqz*M1M2,N1,N2,M1,M2));
      s1z  = real(ffti(dataS1z*M1M2,N1,N2,M1,M2));
      s2z  = real(ffti(dataS2z*M1M2,N1,N2,M1,M2));
      r1z  = real(ffti(dataR1z*M1M2,N1,N2,M1,M2));
      r2z  = real(ffti(dataR2z*M1M2,N1,N2,M1,M2));
      
      %/* Double x-derivatives */
      datahxx=-kx2.*real(datah)-kx2.*imag(datah)*1i;
      dataqxx=-kx2.*real(dataq)-kx2.*imag(dataq)*1i;
      datapxx=-kx2.*real(datap)-kx2.*imag(datap)*1i; 
        
      hxx = real(ffti(datahxx*M1M2,N1,N2,M1,M2));
      pxx = real(ffti(datapxx*M1M2,N1,N2,M1,M2));
      qxx = real(ffti(dataqxx*M1M2,N1,N2,M1,M2));

      
      %/* Double z-derivatives */
      datahzz=-kz2.*real(datah)-kz2.*imag(datah)*1i;
      dataqzz=-kz2.*real(dataq)-kz2.*imag(dataq)*1i;
      datapzz=-kz2.*real(datap)-kz2.*imag(datap)*1i; 
       
      hzz = real(ffti(datahzz*M1M2,N1,N2,M1,M2));
      pzz = real(ffti(datapzz*M1M2,N1,N2,M1,M2));
      qzz = real(ffti(dataqzz*M1M2,N1,N2,M1,M2));

      
      %/* Cross derivatives */
      datahxz=-kxz.*real(datah)-kxz.*imag(datah)*1i;
      dataqxz=-kxz.*real(dataq)-kxz.*imag(dataq)*1i;
      datapxz=-kxz.*real(datap)-kxz.*imag(datap)*1i;

      hxz = real(ffti(datahxz*M1M2,N1,N2,M1,M2));
      pxz = real(ffti(datapxz*M1M2,N1,N2,M1,M2));
      qxz = real(ffti(dataqxz*M1M2,N1,N2,M1,M2));

      
      %/* Gradient of the Laplacien of h (=curvature K) 
      datahLx=-(kx.*k2.*imag(datah)-kx.*k2.*real(datah)*1i);
      datahLz=-(kz.*k2.*imag(datah)-kz.*k2.*real(datah)*1i); 
      
      dhdt=-real(dataqx)-real(datapz)-imag(dataqx)*1i-imag(datapz)*1i;
      
      h   = real(ffti(datah*M1M2,N1,N2,M1,M2));
      P   = real(ffti(datap*M1M2,N1,N2,M1,M2));
      q   = real(ffti(dataq*M1M2,N1,N2,M1,M2));
      s1   = real(ffti(dataS1*M1M2,N1,N2,M1,M2));
      s2   = real(ffti(dataS2*M1M2,N1,N2,M1,M2));
      r1   = real(ffti(dataR1*M1M2,N1,N2,M1,M2));
      r2   = real(ffti(dataR2*M1M2,N1,N2,M1,M2));

datahLx = ffti(datahLx*M1M2,N1,N2,M1,M2);
datahLz = ffti(datahLz*M1M2,N1,N2,M1,M2);

h2=h.*h;
h3=h2.*h;

%The following lines calculate the temporal derivatives for the volume flow
%rate q and p and the corrections s1,s2,r1,r2

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
                18304.*h3.*hZ.*zeta - 18304.*h3.*real(datahLz)))./(delta.*h2) ; ... %CHECK FOR INFLUENCE EL FIELD
                
              
% /* FOURIER SPACE (attention /(M1*M2)!  */
dqdt=fft2(dqdt);
ds1dt=fft2(ds1dt);
ds2dt=fft2(ds2dt);
dpdt=fft2(dpdt);
dr1dt=fft2(dr1dt);
dr2dt=fft2(dr2dt);

%    /*cut to avoid aliasing errors */
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

%Conversion back to the state vector U
dhdt=reshape(dhdt',1,[]);
dqdt=reshape(dqdt',1,[]);
dpdt=reshape(dpdt',1,[]);
ds1dt=reshape(ds1dt',1,[]);
ds2dt=reshape(ds2dt',1,[]);
dr1dt=reshape(dr1dt',1,[]);
dr2dt=reshape(dr2dt',1,[]);

dUdt=[dhdt(1:N/2) dqdt(1:N/2) dpdt(1:N/2) ds1dt(1:N/2) ds2dt(1:N/2) dr1dt(1:N/2) dr2dt(1:N/2)]';
end