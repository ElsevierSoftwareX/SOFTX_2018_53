function [uField,vField,wField,uSurf,vSurf,hx,y,cWave,h,alpha] = reconstruct_velocity_fields(U,M1,M2,N1,N2,h0,kx,kz,ScaleH,ScaleU)

%Calculate the 2d velocity field and the surface velocity for the flow
%Calculate the velocity of the wave for the full model (isothermal flow)
%y defines the resolution in crosswise direction
%x defines the resolution in streamwise direction
%z defines the resolution in spanwise direction
sizeU=size(U(1,:));

if sizeU(2) == M2*(N1+1)*7 %Full second-order model
    reconstru=reshape(U(1,:),(N1+1)*M2,7);  
    
    dataS1 = reconstru(:,4);
    dataS2 = reconstru(:,5);
    dataR1 = reconstru(:,6);
    dataR2 = reconstru(:,7);
else %Simplified model
    reconstru=reshape(U(1,:),(N1+1)*M2,3);  
    
    dataS1 = reconstru(:,1)*0; %Set s1,s2,r1,r2 to zero for the simplified model
    dataS2 = dataS1;
    dataR1 = dataS2;
    dataR2 = dataS2;
end

%reconstru=reshape(yU,(N1+1)*M2,3);
datah = reconstru(:,1);
dataq = reconstru(:,2);
datap = reconstru(:,3);

clear reconstru

datah=reshape(datah,M2,N1+1)';
dataq=reshape(dataq,M2,N1+1)';
datap=reshape(datap,M2,N1+1)';
dataS1=reshape(dataS1,M2,N1+1)';
dataS2=reshape(dataS2,M2,N1+1)';
dataR1=reshape(dataR1,M2,N1+1)';
dataR2=reshape(dataR2,M2,N1+1)';


datah(1,1)=h0+imag(datah(1,1));

% %kx and kz will be defined in the main run
% 
kx2=kx.*kx;
kz2=kz.*kz;
kxz=kx.*kz;
k2=kx2+kz2;
% 
 M1M2=M1*M2;

      %/* x-derivatives */
      datahx=kx.*imag(datah)-kx.*real(datah)*i;      
      dataqx=kx.*imag(dataq)-kx.*real(dataq)*i;
      datapx=kx.*imag(datap)-kx.*real(datap)*i;
      dataS1x=kx.*imag(dataS1)-kx.*real(dataS1)*i;
      dataS2x=kx.*imag(dataS2)-kx.*real(dataS2)*i;
      dataR1x=kx.*imag(dataR1)-kx.*real(dataR1)*i;
      dataR2x=kx.*imag(dataR2)-kx.*real(dataR2)*i;
      
      hx  = real(ffti(datahx*M1M2,N1,N2,M1,M2));
      px  = real(ffti(datapx*M1M2,N1,N2,M1,M2));
      qx  = real(ffti(dataqx*M1M2,N1,N2,M1,M2));
      s1x  = real(ffti(dataS1x*M1M2,N1,N2,M1,M2));
      s2x  = real(ffti(dataS2x*M1M2,N1,N2,M1,M2));
      r1x  = real(ffti(dataR1x*M1M2,N1,N2,M1,M2));
      r2x  = real(ffti(dataR2x*M1M2,N1,N2,M1,M2));
      
      
      %/* z-derivatives */
      datahz=kz.*imag(datah)-kz.*real(datah)*i;
      dataqz=kz.*imag(dataq)-kz.*real(dataq)*i;
      datapz=kz.*imag(datap)-kz.*real(datap)*i;
      dataS1z=kz.*imag(dataS1)-kz.*real(dataS1)*i;
      dataS2z=kz.*imag(dataS2)-kz.*real(dataS2)*i;
      dataR1z=kz.*imag(dataR1)-kz.*real(dataR1)*i;
      dataR2z=kz.*imag(dataR2)-kz.*real(dataR2)*i;
      
      hZ  = real(ffti(datahz*M1M2,N1,N2,M1,M2));
      pz  = real(ffti(datapz*M1M2,N1,N2,M1,M2));
      qz  = real(ffti(dataqz*M1M2,N1,N2,M1,M2));
      s1z  = real(ffti(dataS1z*M1M2,N1,N2,M1,M2));
      s2z  = real(ffti(dataS2z*M1M2,N1,N2,M1,M2));
      r1z  = real(ffti(dataR1z*M1M2,N1,N2,M1,M2));
      r2z  = real(ffti(dataR2z*M1M2,N1,N2,M1,M2));
      
      %/* Double x-derivatives */
      datahxx=-kx2.*real(datah)-kx2.*imag(datah)*i;
      dataqxx=-kx2.*real(dataq)-kx2.*imag(dataq)*i;
      datapxx=-kx2.*real(datap)-kx2.*imag(datap)*i; 
      
      
      hxx = real(ffti(datahxx*M1M2,N1,N2,M1,M2));
      pxx = real(ffti(datapxx*M1M2,N1,N2,M1,M2));
      qxx = real(ffti(dataqxx*M1M2,N1,N2,M1,M2));

      
      %/* Double z-derivatives */
      datahzz=-kz2.*real(datah)-kz2.*imag(datah)*i;
      dataqzz=-kz2.*real(dataq)-kz2.*imag(dataq)*i;
      datapzz=-kz2.*real(datap)-kz2.*imag(datap)*i; 

          
      hzz = real(ffti(datahzz*M1M2,N1,N2,M1,M2));
      pzz = real(ffti(datapzz*M1M2,N1,N2,M1,M2));
      qzz = real(ffti(dataqzz*M1M2,N1,N2,M1,M2));

      
      %/* Cross derivatives */  %CHECK VORZEICHEN
      datahxz=-kxz.*real(datah)-kxz.*imag(datah)*i;
      dataqxz=-kxz.*real(dataq)-kxz.*imag(dataq)*i;
      datapxz=-kxz.*real(datap)-kxz.*imag(datap)*i;

      
      hxz = real(ffti(datahxz*M1M2,N1,N2,M1,M2));
      pxz = real(ffti(datapxz*M1M2,N1,N2,M1,M2));
      qxz = real(ffti(dataqxz*M1M2,N1,N2,M1,M2));

      
      %/* Gradient of the laplacien of h (=curvature K) 
      datahLx=-(kx.*k2.*imag(datah)-kx.*k2.*real(datah)*i);
      datahLz=-(kz.*k2.*imag(datah)-kz.*k2.*real(datah)*i); 
      
      dhdt=-real(dataqx)-real(datapz)-imag(dataqx)*i-imag(datapz)*i;
      
h   = real(ffti(datah*M1M2,N1,N2,M1,M2));
P   = real(ffti(datap*M1M2,N1,N2,M1,M2));
q   = real(ffti(dataq*M1M2,N1,N2,M1,M2));
s1   = real(ffti(dataS1*M1M2,N1,N2,M1,M2));
s2   = real(ffti(dataS2*M1M2,N1,N2,M1,M2));
r1   = real(ffti(dataR1*M1M2,N1,N2,M1,M2));
r2   = real(ffti(dataR2*M1M2,N1,N2,M1,M2));

hmax=max(max(h));
 y=0:0.01*hmax:hmax*1.2;   
  y=y*ScaleH;            

for j=1:1:M1
    for k=1:1:M2
%Determine yBar with datah in real space
yBar(j,k,:)=(y./h(j,k));
    end
end

%Calculate functions F0-F3 and G0-G3

F0=yBar-1/2*yBar.^2;

F1=yBar-17/6.*yBar.^2+7/3.*yBar.^3-7/12.*yBar.^4;

F2=yBar-13/2.*yBar.^2+57/4*yBar.^3-111/8.*yBar.^4+99/16.*yBar.^5-33/32.*yBar.^6;

F3=yBar-531/62.*yBar.^2+2871/124*yBar.^3-6369/248.*yBar.^4+29601/2480.*yBar.^5-9867/4960.*yBar.^6;


[k l]=size(q);
for k=1:1:l
    F0Star=squeeze(F0(:,k,:));
    F1Star=squeeze(F1(:,k,:));
    F2Star=squeeze(F2(:,k,:));
    [j m]=size(y);
    for j=1:1:m
        uField(:,k,j)=3*(q(:,k)-s1(:,k)-s2(:,k))./h(:,k).*F0Star(:,j)+45*s1(:,k)./h(:,k).*F1Star(:,j)+210*s2(:,k)./h(:,k).*F2Star(:,j);
        wField(:,k,j)=3*(P(:,k)-r1(:,k)-r2(:,k))./h(:,k).*F0Star(:,j)+45*r1(:,k)./h(:,k).*F1Star(:,j)+210*r2(:,k)./h(:,k).*F2Star(:,j);
    
        
     
       vField(:,k,j)=(y(:,j).^2.*(-24.*h(:,k).^6.*(pz(:,k) + qx(:,k) + 14.*r1z(:,k) + 69.*r2z(:,k) + 14.*s1x(:,k) + 69.*s2x(:,k)) - 3465.*(hZ(:,k).*r2(:,k) + hx(:,k).*s2(:,k)).*y(:,j).^5 + 495.*h(:,k).*y(:,j).^4.*(42.*hZ(:,k).*r2(:,k) + 42.*hx(:,k).*s2(:,k) + r2z(:,k).*y(:,j) + s2x(:,k).*y(:,j)) - ...
          105.*h(:,k).^2.*y(:,j).^3.*(4.*hZ(:,k).*r1(:,k) + 444.*hZ(:,k).*r2(:,k) + 4.*hx(:,k).*s1(:,k) + 444.*hx(:,k).*s2(:,k) + 33.*r2z(:,k).*y(:,j) + 33.*s2x(:,k).*y(:,j)) + ...
          84.*h(:,k).^3.*y(:,j).^2.*(20.*hZ(:,k).*r1(:,k) + 570.*hZ(:,k).*r2(:,k) + 20.*hx(:,k).*s1(:,k) + 570.*hx(:,k).*s2(:,k) + r1z(:,k).*y(:,j) + 111.*r2z(:,k).*y(:,j) + s1x(:,k).*y(:,j) + 111.*s2x(:,k).*y(:,j)) + ...
          8.*h(:,k).^5.*(6.*hZ(:,k).*P(:,k) + 6.*hx(:,k).*q(:,k) + 84.*hZ(:,k).*r1(:,k) + 414.*hZ(:,k).*r2(:,k) + 84.*hx(:,k).*s1(:,k) + 414.*hx(:,k).*s2(:,k) + pz(:,k).*y(:,j) + qx(:,k).*y(:,j) + 84.*r1z(:,k).*y(:,j) + 909.*r2z(:,k).*y(:,j) + 84.*s1x(:,k).*y(:,j) + 909.*s2x(:,k).*y(:,j)) - ...
          6.*h(:,k).^4.*y(:,j).*(4.*hZ(:,k).*P(:,k) + 4.*hx(:,k).*q(:,k) + 336.*hZ(:,k).*r1(:,k) + 3636.*hZ(:,k).*r2(:,k) + 336.*hx(:,k).*s1(:,k) + 3636.*hx(:,k).*s2(:,k) + 70.*r1z(:,k).*y(:,j) + 1995.*r2z(:,k).*y(:,j) + 70.*s1x(:,k).*y(:,j) + 1995.*s2x(:,k).*y(:,j))))./(16..*h(:,k).^8);         
    
    
    end
end


uField=uField*ScaleU;
vField=vField*ScaleU;
wField=wField*ScaleU;
h=h*ScaleH;

[k l m]=size(uField);

for j=1:1:M1
    for k=1:1:M2
        for l=1:1:m
             if y(l)>h(j,k)
                 uField(j,k,l)=nan;
                 vField(j,k,l)=nan;
                 wField(j,k,l)=nan;
                 alpha(j,k,l)=0;
             else
                 alpha(j,k,l)=1;
             end
        end
    end
end


yBar=1;

F0=yBar-1/2*yBar.^2;

F1=yBar-17/6.*yBar.^2+7/3.*yBar.^3-7/12.*yBar.^4;

F2=yBar-13/2.*yBar.^2+57/4*yBar.^3-111/8.*yBar.^4+99/16.*yBar.^5-33/32.*yBar.^6;

F3=yBar-531/62.*yBar.^2+2871/124*yBar.^3-6369/248.*yBar.^4+29601/2480.*yBar.^5-9867/4960.*yBar.^6;

[ac l]=size(q);

for k=1:1:l
    F0Star=F0;
    F1Star=F1;
    F2Star=F2;
    
    [j m]=size(y);
    
    j=1;
       uSurf(:,k)=3*(q(:,k)-s1(:,k)-s2(:,k))./h(:,k).*F0Star(:,j)+45*s1(:,k)./h(:,k).*F1Star(:,j)+210*s2(:,k)./h(:,k).*F2Star(:,j);
       wSurf(:,k)=3*(P(:,k)-r1(:,k)-r2(:,k))./h(:,k).*F0Star(:,j)+45*r1(:,k)./h(:,k).*F1Star(:,j)+210*r2(:,k)./h(:,k).*F2Star(:,j);

       vSurf(:,k)=(h(:,k).^2.*(-24.*h(:,k).^6.*(pz(:,k) + qx(:,k) + 14.*r1z(:,k) + 69.*r2z(:,k) + 14.*s1x(:,k) + 69.*s2x(:,k)) - 3465.*(hZ(:,k).*r2(:,k) + hx(:,k).*s2(:,k)).*h(:,k).^5 + 495.*h(:,k).*h(:,k).^4.*(42.*hZ(:,k).*r2(:,k) + 42.*hx(:,k).*s2(:,k) + r2z(:,k).*h(:,k) + s2x(:,k).*h(:,k)) - ...
          105.*h(:,k).^2.*h(:,k).^3.*(4.*hZ(:,k).*r1(:,k) + 444.*hZ(:,k).*r2(:,k) + 4.*hx(:,k).*s1(:,k) + 444.*hx(:,k).*s2(:,k) + 33.*r2z(:,k).*h(:,k) + 33.*s2x(:,k).*h(:,k)) + ...
          84.*h(:,k).^3.*h(:,k).^2.*(20.*hZ(:,k).*r1(:,k) + 570.*hZ(:,k).*r2(:,k) + 20.*hx(:,k).*s1(:,k) + 570.*hx(:,k).*s2(:,k) + r1z(:,k).*h(:,k) + 111.*r2z(:,k).*h(:,k) + s1x(:,k).*h(:,k) + 111.*s2x(:,k).*h(:,k)) + ...
          8.*h(:,k).^5.*(6.*hZ(:,k).*P(:,k) + 6.*hx(:,k).*q(:,k) + 84.*hZ(:,k).*r1(:,k) + 414.*hZ(:,k).*r2(:,k) + 84.*hx(:,k).*s1(:,k) + 414.*hx(:,k).*s2(:,k) + pz(:,k).*h(:,k) + qx(:,k).*h(:,k) + 84.*r1z(:,k).*h(:,k) + 909.*r2z(:,k).*h(:,k) + 84.*s1x(:,k).*h(:,k) + 909.*s2x(:,k).*h(:,k)) - ...
          6.*h(:,k).^4.*h(:,k).*(4.*hZ(:,k).*P(:,k) + 4.*hx(:,k).*q(:,k) + 336.*hZ(:,k).*r1(:,k) + 3636.*hZ(:,k).*r2(:,k) + 336.*hx(:,k).*s1(:,k) + 3636.*hx(:,k).*s2(:,k) + 70.*r1z(:,k).*h(:,k) + 1995.*r2z(:,k).*h(:,k) + 70.*s1x(:,k).*h(:,k) + 1995.*s2x(:,k).*h(:,k))))./(16..*h(:,k).^8); 
      for ls = 1:ac
          if (hx(ls,k)>0.005)
           cWave(ls,k)=(hx(ls,k).*uSurf(ls,k)-vSurf(ls,k))./hx(ls,k);
          else
           cWave(ls,k)=NaN;
          end
     end
end

end
