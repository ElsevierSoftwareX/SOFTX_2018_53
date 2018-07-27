% *******************************************************************
% *                        WAVEMAKER                                *
% * Software to simulate falling films using the Integral boundary
% * layer method ....
% *   
% *                                                                 *
% *                   Wilko Rohlfs, Manuel Rietz                    *
% *  Institute of Heat and Mass Transfer, RWTH Aachen University    *
% *                       and Benoit Scheid
% *       TIPs lab – Microfluidics, Université Libre de Bruxelles   *
% *******************************************************************
%
% Copyright 2017 Wilko Rohlfs, Manuel Rietz and Benoit Scheid
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% This script is intended to conduct the simulations by
%   Choosing GPU mode or not (if GPU is available)
%   

function  [U,ttime] = runSimulation(handles,hObjectToggle,parameter)

format long
        warning('off');
        mkdir(parameter.Name);
        warning('on');
        save(['./' parameter.Name '/ParameterFile'],'parameter');
       
        %Choose the fft method. Used to optimize fft speed
        fftw('planner','patient');
        
        %Check if the machine supports GPU.
        %GPU is used if the machine supports it.
        try
            d = gpuDevice;
            if parameter.GPU_Acceleration == true
            GPU = true;
            else
            GPU = false;    
            end
        catch
            GPU = false;
        end

        %Accuracy of the simulation Fourier coefficients in streamwise and spanwiese direction
        NN1=parameter.NNst  ;
        NN2=parameter.NNsp  ;
               
        %Read dimensionless parameters in Shkadov scaling
        delta = parameter.Shkadov_delta;
        zeta = parameter.Shkadov_zeta;
        eta = parameter.Shkadov_eta;
        
        k0x = parameter.Shkadov_k0x;
        k0z = parameter.Shkadov_k0z;
        h0=parameter.h0;
         
        %Read maximum time and print time step
        tmax = str2num(handles.TMax.String);
        prtstep = str2num(handles.PrtStep.String);
        
        %Determine viscous time scale tv and viscous length scale lv and Nusselt
        %film thickness hn
        DimVar_gx=parameter.DimVar_gAbs*sin(parameter.DimVar_Theta/360*2*pi());
        lv=(parameter.DimVar_nu^2/DimVar_gx)^(1/3);
        tv=(parameter.DimVar_nu/(DimVar_gx^2))^(1/3);
        hnDach=(3*parameter.Reynolds)^(1/3)*lv;
        
        %Determine tmax and PrintTimestep for different scalings depending
        %on user choice
        if handles.New_LstBxValue == 1
           tmax=tmax*hnDach/(tv*lv*parameter.Shkadov_kappa);
           prtstep=prtstep*hnDach/(tv*lv*parameter.Shkadov_kappa);
        end
        
        if handles.New_LstBxValue == 2
           tmax=tmax/parameter.Shkadov_kappa;
           prtstep=prtstep/parameter.Shkadov_kappa;
        end
        
             
        N1=2^NN1; M1=2*N1; N2=2^NN2; M2=2*N2; N=4*N2*(N1+1);

        %kappa=1/sqrt(eta);
           
        if(k0x~=0) Lx=2.*pi/k0x; end
        if(k0z~=0) Lz=2.*pi/k0z; end
        if(NN1==0) k0x=0;Lx=0; end
        if(NN2==0) k0z=0;Lz=0; end
        if(NN1==0||NN2==0) dim=2; else dim=3; end
        
        AlaisingRadius = parameter.AlaisingRadius;  %Radius for anti-alaising
        %Create Alaising filter and Initial conditions for CPU
        [kx,kz,k2cut]=CPU_AliasingFilter(AlaisingRadius,N1,N2,M1,M2,k0x,k0z);

        if parameter.continue==0
        %U is the space variable which containes h,q,...
        U = Create_Initial_conditions(parameter,h0,Lz,Lx,k0x,k0z,NN1,NN2,k2cut,N);
        ttime = 0; %Start with time 0
        else
            U = parameter.UEnd;
            ttime = parameter.TEnd;
        end
        
        %Create Alaising filter for GPU
        if(GPU == true)
                [kxGPU,kzGPU,k2cut,NoConj,Conj]=GPU_AliasingFilter(AlaisingRadius,N1,N2,M1,M2,k0x,k0z);
        end
        
        %Initialize variables for storage and evaluation
        storeData.time = [];
        storeData.hmax = [];
        storeData.hmin = [];
        storeData.q = [];
               
        while(tmax-ttime>prtstep && get(hObjectToggle,'Value')==1)
            
            stop = ttime+prtstep;
            while (ttime<stop)        % Loop over prtstep
                
                % Fifth-order Runge-Kutta with adaptative stepsize */
                options = odeset('RelTol',1e-3);
                tic
                T = [ttime ttime+prtstep/2 ttime+prtstep]; %Initiate three steps otherwise ode45 returns all intermediate time step results
                if parameter.Model == 1 %Start with simplified model
                    if(GPU == true)
                                             
                        [T,U] = ode45(@(T,U)GPU_dUdt_SimplifiedModel(U,h0,M1,M2,N1,N2,kxGPU,kzGPU,delta,eta,zeta,k2cut,N,Conj,NoConj,GPU),T,U,options);    
                        
                    else
                        [T,U] = ode45(@(T,U)CPU_dUdt_SimplifiedModel(U,h0,M1,M2,N1,N2,kx,kz,delta,eta,zeta,k2cut,N,0,0),T,U,options);        
                    end
                elseif parameter.Model == 2 %Start Full second-order model
                    if(GPU == true)
                        [T,U] = ode45(@(T,U)GPU_dUdt_FullModel(U,h0,M1,M2,N1,N2,kxGPU,kzGPU,delta,eta,zeta,k2cut,N,Conj,NoConj,GPU),T,U,options);
                    else
                        [T,U] = ode45(@(T,U)CPU_dUdt_FullModel(U,h0,M1,M2,N1,N2,kx,kz,delta,eta,zeta,k2cut,N),T,U,options);
                    end
                end
                
                clear temp
                ttime = T(end);
                U = U(end,:);
                [storeData] = users_output(handles,parameter,U,M1,M2,N1,N2,NN2,h0,storeData,kx,kz,ttime);
        tout=ttime;        
        if handles.New_LstBxValue == 1
           tout=tout/hnDach*(tv*lv*parameter.Shkadov_kappa);
        end
        
        if handles.New_LstBxValue == 2
           tout=tout*parameter.Shkadov_kappa;
        end
        
                disp(['Current timestep: ' num2str(tout)])
                toc
            end
        end


    end