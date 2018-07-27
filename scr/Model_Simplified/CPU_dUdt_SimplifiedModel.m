function [dUdt] = CPU_dUdt_SimplifiedModel(yU,h0,M1,M2,N1,N2,kx,kz,delta,eta,zeta,k2cut,N,A,hs)
%function calculating the equation for the simplified model that will be solved by ode45
%
%see also fft, ode45
% 

%reconstruct data from U-vector (filmthickness h, streamwise flow rate q
%and spanwise flow rate p

reconstru=reshape(yU,(N1+1)*M2,3);
datah = reconstru(:,1);
dataq = reconstru(:,2);
datap = reconstru(:,3);
clear reconstru

datah=reshape(datah,M2,N1+1)';
dataq=reshape(dataq,M2,N1+1)';
datap=reshape(datap,M2,N1+1)';

%fix the value of the average film thickness
datah(1,1)=h0+imag(datah(1,1)); %Set average film thickness to fixed value

%Wavenumber values needed for derivatives in Fourier space
kx2=kx.*kx;
kz2=kz.*kz;
kxz=kx.*kz;
k2=kx2+kz2;

M1M2=M1*M2;

      %/* x-derivatives */
datahx=kx.*imag(datah)-kx.*real(datah)*1i;      
dataqx=kx.*imag(dataq)-kx.*real(dataq)*1i;
datapx=kx.*imag(datap)-kx.*real(datap)*1i;

%Transform variables back to real space
hx  = real(ffti(datahx*M1M2,N1,N2,M1,M2));
px  = real(ffti(datapx*M1M2,N1,N2,M1,M2));
qx  = real(ffti(dataqx*M1M2,N1,N2,M1,M2));

      
      %/* z-derivatives */
datahz=kz.*imag(datah)-kz.*real(datah)*1i;
dataqz=kz.*imag(dataq)-kz.*real(dataq)*1i;
datapz=kz.*imag(datap)-kz.*real(datap)*1i;

%Transform variables back to real space
hZ  = real(ffti(datahz*M1M2,N1,N2,M1,M2));
pz  = real(ffti(datapz*M1M2,N1,N2,M1,M2));
qz  = real(ffti(dataqz*M1M2,N1,N2,M1,M2));

      
      %/* Double x-derivatives */
datahxx=-kx2.*real(datah)-kx2.*imag(datah)*1i;
dataqxx=-kx2.*real(dataq)-kx2.*imag(dataq)*1i;
datapxx=-kx2.*real(datap)-kx2.*imag(datap)*1i;  

%Transform variables back to real space      
hxx = real(ffti(datahxx*M1M2,N1,N2,M1,M2));
pxx = real(ffti(datapxx*M1M2,N1,N2,M1,M2));
qxx = real(ffti(dataqxx*M1M2,N1,N2,M1,M2));
   
      %/* Double z-derivatives */
datahzz=-kz2.*real(datah)-kz2.*imag(datah)*1i;
dataqzz=-kz2.*real(dataq)-kz2.*imag(dataq)*1i;
datapzz=-kz2.*real(datap)-kz2.*imag(datap)*1i; 

%Transform variables back to real space
hzz = real(ffti(datahzz*M1M2,N1,N2,M1,M2));
pzz = real(ffti(datapzz*M1M2,N1,N2,M1,M2));
qzz = real(ffti(dataqzz*M1M2,N1,N2,M1,M2));

      
     %/* Cross derivatives */
datahxz=-kxz.*real(datah)-kxz.*imag(datah)*1i;
dataqxz=-kxz.*real(dataq)-kxz.*imag(dataq)*1i;
datapxz=-kxz.*real(datap)-kxz.*imag(datap)*1i;

%Transform variables back to real space
hxz = real(ffti(datahxz*M1M2,N1,N2,M1,M2));
pxz = real(ffti(datapxz*M1M2,N1,N2,M1,M2));
qxz = real(ffti(dataqxz*M1M2,N1,N2,M1,M2));

    %/* Gradient of the laplacien of h (=curvature K) 
datahLx=-(kx.*k2.*imag(datah)-kx.*k2.*real(datah)*1i);
datahLz=-(kz.*k2.*imag(datah)-kz.*k2.*real(datah)*1i); 
    
dhdt=-real(dataqx)-real(datapz)-imag(dataqx)*1i-imag(datapz)*1i;

%Transform variables back to real space
h   = real(ffti(datah*M1M2,N1,N2,M1,M2));
P   = real(ffti(datap*M1M2,N1,N2,M1,M2));
q   = real(ffti(dataq*M1M2,N1,N2,M1,M2));


datahLx = ffti(datahLx*M1M2,N1,N2,M1,M2);
datahLz = ffti(datahLz*M1M2,N1,N2,M1,M2);

h2=h.*h;
h3=h2.*h;

Ineq=0; %Switch to 1 for second order corrections
       
dqdt = (0.011904761904761904.*(8.*h.*(4900.*A.*hs.*hx - 735.*eta.*(13.*hx.*hZ.*P + 16.*hx.*hx.*q + 3.*hZ.*hZ.*q) + ...
                  6.*q.*(1225. + 9.*delta.*(-70.*hZ.*P + hx.*q.*(-70. + delta.*hx.*Ineq.*q)))) + ...
                  6.*h2.*(245.*eta.*(73.*hxz.*P + 43.*hZ.*px + 13.*hx.*pz + 96.*hxx.*q + 23.*hzz.*q + 72.*hx.*qx + ...
                  16.*hZ.*qz) + 8.*delta.*(560.*pz.*q - 17.*q.*(-70. + delta.*hx.*Ineq.*q).*qx + 630.*P.*qz)) + ...
                  - 5880.*h3.*(2.*eta.*(7.*pxz + 9.*qxx + 2.*qzz)) - ...
                  300.*hx.*(196.*A.*hs.*hs) - ...
                  35.*h2.*h2.*(560. - 560.*hx.*zeta) - ...
                  19600.*h2.*h2.*real(datahLx)))./(delta.*h2.*(4.*h.*(-70. + delta.*hx.*Ineq.*q)));

                  
dpdt = (0.002976190476190476.*(-280.*A.*(2.*h - 3.*hs).*hs.*hZ + ...
                   h.*(12.*(P.*(-70. + 7.*eta.*(3.*hx.*hx + 16.*hZ.*hZ) + 36.*delta.*hZ.*P) + ...
                   hx.*(91.*eta.*hZ + 36.*delta.*P).*q) - 3.*h.*(16.*delta.*(17.*P.*pz + 9.*px.*q + 8.*P.*qx) + ...
                   7.*eta.*(23.*hxx.*P + 96.*hzz.*P + 16.*hx.*px + 72.*hZ.*pz + 73.*hxz.*q + 13.*hZ.*qx + 43.*hx.*qz)) + ...
                   84.*h2.*(2.*eta.*(2.*pxx + 9.*pzz + 7.*qxz)) - 280.*h3.*hZ.*zeta ) + ...
                   280.*h2.*h2.*real(datahLz)))./(delta.*h3);
                          
% /* FOURIER SPACE (attention /(M1*M2)!  */
dqdt=fft2(dqdt);
dpdt=fft2(dpdt);

%    /* Alaising filter removes high frequency modes */
dqdt=dqdt.*k2cut/M1M2;
dpdt=dpdt.*k2cut/M1M2;

%Transform data back to state vector
dhdt=conj(dhdt);
dqdt=conj(dqdt);
dpdt=conj(dpdt);

dhdt=reshape(dhdt',1,[]);
dqdt=reshape(dqdt',1,[]);
dpdt=reshape(dpdt',1,[]);

%Required derivative 
dUdt=[dhdt(1:N/2) dqdt(1:N/2) dpdt(1:N/2)]';

end