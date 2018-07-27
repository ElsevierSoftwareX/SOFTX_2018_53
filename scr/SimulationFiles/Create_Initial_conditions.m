%This function creates the initial condition depending on the given
%parameter
%Lx and Lz:
%k0x and k0z:
%NN1 and NN2:
%k2cut:
%N:
%Ax Bx Az Bz:
%nwx and nwz:

function [U] = Create_Initial_conditions(parameter,h0,Lz,Lx,k0x,k0z,NN1,NN2,k2cut,N)

if parameter.UseExistingInitialCondition==false
    Ax=parameter.stSine; Az=parameter.spSine; Bx=parameter.stCosine; Bz=parameter.spCosine; nwx=parameter.stNwaves; nwz=parameter.spNwaves; ampl=parameter.stNoise; amplz=parameter.spNoise;
    
    if(NN1==0) k0x=0;Lx=0;Ax=0;Bx=0;NN1P=1; end
    if(NN2==0) k0z=0;Lz=0;Az=0;Bz=0;NN2P=1; end
    
    [datah,datap,dataq,U] = Initial_condition(ampl, Lz, Lx, h0, Ax, k0x, Bx, nwx, k0z, Az, nwz, Bz, amplz, NN1, NN2);
    
    datas1s2r1r2=reshape(datah',1,[])*0;
    
    %Remove high frequency parts by aliasing filter
    datah=datah.*k2cut;
    datap=datap.*k2cut;
    dataq=dataq.*k2cut;
    
    %Reshape data for U vector
    datahSV=reshape(datah',1,[]);
    datapSV=reshape(datap',1,[]);
    dataqSV=reshape(dataq',1,[]);
    
    U=[datahSV(1:N/2) dataqSV(1:N/2) datapSV(1:N/2)];
    clear datahSV datapSV dataqSV dataTSV
       
    if parameter.Model==2 %For full second order model add dummy values for s1,s2,r1,r2
    U=[U datas1s2r1r2(1:N/2) datas1s2r1r2(1:N/2) datas1s2r1r2(1:N/2) datas1s2r1r2(1:N/2)]; %Set initial conditions for s1 s2 r1 r2
    end
    clear amplz Ax Bx Az Bz nwz nwx ampl datas1s2r1r2
end

%Load initial condition
if parameter.UseExistingInitialCondition==true
    fileName=[parameter.LoadInitialConditionPathName '/' parameter.LoadInitialConditionFileName];
    U=load(fileName);
    U=U.Usave;
end
end