function [storeData] = users_output(handles,parameter,U,M1,M2,N1,N2,NN2,h0,storeData,kx,kz,ttime)
%make all the output for the user

[datah,dataq,~]   = reconstruct_U(U,N1,N2,M1,M2);
hmax=real(max(max(datah)));
hmin=real(min(min(datah)));
q = mean(mean(real(dataq)*parameter.Reynolds*3));
storeData.hmax=[storeData.hmax hmax];
storeData.hmin=[storeData.hmin hmin];
storeData.time=[storeData.time ttime];
storeData.q=[storeData.q q];

%Introduce Scalings for plots

if parameter.Scaling==1 %Dimensional scaling
    ScaleX=parameter.DimVar_Lx;
    ScaleZ=parameter.DimVar_Lz;
    DimVar_gx=parameter.DimVar_gAbs*sin(parameter.DimVar_Theta/360*2*pi());
    lv=(parameter.DimVar_nu^2/DimVar_gx)^(1/3);
    tv=(parameter.DimVar_nu/(DimVar_gx^2))^(1/3);
    hnDach=(3*parameter.Reynolds)^(1/3)*lv;
    ScaleH=hnDach;
    ScaleQ=1;  % The volume flow rate is still described by Re for convinience
    ScaleT=tv*lv/hnDach*parameter.Shkadov_kappa;
    ScaleU=hnDach^2/lv/tv;
elseif parameter.Scaling==2 %Nusselt scaling
    ScaleX = 2*pi()/parameter.Nusselt_k0x;
    ScaleZ = 2*pi()/parameter.Nusselt_k0z;
    ScaleH = 1;
    ScaleQ = 1;
    ScaleT=parameter.Shkadov_kappa;
    ScaleU = 1;
elseif parameter.Scaling==3 %Shkadov scaling
    ScaleX = 2*pi()/parameter.Shkadov_k0x;
    ScaleZ = 2*pi()/parameter.Shkadov_k0z;
    ScaleH = 1;
    ScaleQ = parameter.Shkadov_delta/parameter.Reynolds;
    ScaleT = 1;
    ScaleU = 1;
end

set(handles.Time,'String',['Time: ' num2str(ttime*ScaleT)]);

if(NN2==0) %Plot for 2-D cases
    
cla(handles.axes2,'reset')
    
    hold(handles.axes2,'on');
    axis(handles.axes2,[0 1 0 max(storeData.hmax)*1.1]);
    [vU,~,~,~,~,~,y,cWave,~,~]= reconstruct_velocity_fields(U(end,:),M1,M2,N1,N2,h0,kx,kz,1,1);
    vU=squeeze(vU(:,1,:))';

    if (nanstd(cWave)<0.1)
        psi=cumsum(vU-nanmedian(cWave(:,1)));
        set(handles.Result_CWave,'String',['Wave velocity: ' num2str(nanmedian(cWave(:,1)*ScaleU))]);
        set(handles.Result_fWave,'String',['Wave frequency: ' num2str(nanmedian(cWave(:,1)*ScaleU)/parameter.DimVar_Lx)]);
        set(handles.FrameOfReference,'String','Moving frame of reference');
    else
        psi=cumsum(vU);
        set(handles.Result_CWave,'String','Wave velocity: non stationary');
        set(handles.Result_fWave,'String','Wave frequency: non stationary');
        set(handles.FrameOfReference,'String','Stationary frame of reference');
    end
    contour(handles.axes2,(0:1/(M1-1):1)*ScaleX,y*ScaleH,fliplr(-psi),20);
    contour(handles.axes2,(0:1/(M1-1):1)*ScaleX,y*ScaleH,-psi,[0 1]);
    plot(handles.axes2,(0:1/(M1-1):1)*ScaleX,flipud(real(datah(:,1)))*ScaleH,'r','LineWidth',2);
    axis(handles.axes2,[0 ScaleX 0 max(storeData.hmax)*ScaleH*1.1]);
    if parameter.Scaling==1
        xlabel(handles.axes2,'Streamwise coordinate [m]')
        ylabel(handles.axes2,'Film thickness [m]')
    else
        ylabel(handles.axes2,'Film thickness')
        xlabel(handles.axes2,'Streamwise coordinate ')
    end
    box(handles.axes2,'on')
    drawnow
    
else %Plot for three dimensional cases

    hold(handles.axes2,'off');
    set(handles.Result_CWave,'String','');
    set(handles.Result_fWave,'String','');
    set(handles.FrameOfReference,'String','');
    surf(handles.axes2,(0:1/(M2-1):1)*ScaleZ,(0:1/(M1-1):1)*ScaleX,real(datah)*ScaleH,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    view(handles.axes2,[parameter.azFig1 parameter.elFig1]);
    %camproj(handles.axes2,'perspective');
    axis(handles.axes2,[0 ScaleZ 0 ScaleX min(storeData.hmin)*ScaleH max(storeData.hmax)*ScaleH]);

end

%Plot in the lower right figure

% check for matlab version regarding use of yyaxis or plotyy (R2015b and older)
if ~verLessThan('matlab', '9.0')            

view(handles.axes3,2)
yyaxis(handles.axes3,'left')

plot(handles.axes3,storeData.time*ScaleT,storeData.hmax*ScaleH,'b-','MarkerFaceColor','b')
hold(handles.axes3,'on')
plot(handles.axes3,storeData.time*ScaleT,storeData.hmin*ScaleH,'b--')

if parameter.Scaling==1
    ylabel(handles.axes3,'Film thickness [m]')
    xlabel(handles.axes3,'Time [s]')
else
    ylabel(handles.axes3,'Film thickness')
    xlabel(handles.axes3,'Time')
end

yyaxis(handles.axes3,'right')
hold(handles.axes3,'off')
plot(handles.axes3,storeData.time*ScaleT,storeData.q*ScaleQ,'r-','MarkerFaceColor','r')
% view(handles.axes3,[min(storeData.time) max(storeData.time)]);

if parameter.Scaling==3
    ylabel(handles.axes3,'Reduced Reynolds number, \delta_q')   
else
    ylabel(handles.axes3,'Reynolds number, Re_q')
end


drawnow

else

view(handles.axes3,2)

ax=plotyy(handles.axes3,[storeData.time*ScaleT,storeData.time*ScaleT],[storeData.hmax*ScaleH,storeData.hmin*ScaleH],storeData.time*ScaleT,storeData.q*ScaleQ);

if parameter.Scaling==1
    ylabel(ax(1),'Film thickness [m]')
    xlabel(ax(1),'Time [s]')
else
    ylabel(ax(1),'Film thickness')
    xlabel(ax(1),'Time')
end

if parameter.Scaling==3
    ylabel(ax(2),'Reduced Reynolds number, \delta_q')   
else
    ylabel(ax(2),'Reynolds number, Re_q')
end  

drawnow
end
%% save data to specified data formats

% get chosen data formats to save to
handles.saveoptions=getappdata(0,'saveoptions');
[vU,vV,vW,~,~,~,~,~,h,alpha] = reconstruct_velocity_fields(U,M1,M2,N1,N2,h0,kx,kz,ScaleH,ScaleU);

if handles.saveoptions.vtk ==1
    
    %% create 3D dataset
    
    hmax=max(max(h));
    z1=0:0.01*hmax:hmax*1.2;
    x1=0:parameter.DimVar_Lx/(M1-1):parameter.DimVar_Lx;       % streamwise extend of domain
    y1=0:parameter.DimVar_Lz/(M2-1):parameter.DimVar_Lz;       % spanwise extend of domain
    
    [x, y, z] = meshgrid(x1,y1,z1);
    
    u=vU;
    v=vV;
    w=vW;
    
    %% write .vtk file
    
    %Output file name
    filename=([parameter.Name '/' sprintf('%06d',round(ttime*ScaleT,3)) 'Velocity_field.vtk']);
    
    nr_of_elements=numel(x);
    fid = fopen(filename, 'w');
    
    %ASCII file header
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'VTK from Matlab\n');
    fprintf(fid, 'BINARY\n\n');
    fprintf(fid, 'DATASET STRUCTURED_GRID\n');
    fprintf(fid, ['DIMENSIONS ' num2str(size(x,1)) ' ' num2str(size(x,2)) ' ' num2str(size(x,3)) '\n']);
    fprintf(fid, ['POINTS ' num2str(nr_of_elements) ' float\n']);
    fclose(fid);
    
    %append binary x,y,z data
    fid = fopen(filename, 'a');
    fwrite(fid, [reshape(x,1,nr_of_elements);  reshape(y,1,nr_of_elements); reshape(z,1,nr_of_elements)],'float','b');
    
    %append another ASCII sub header
    fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements) '\n']);
    fprintf(fid, 'VECTORS velocity_vectors float\n');
    
    %append binary u,v,w data
    v_vtk=permute(v,[2 1 3]);
    u_vtk=permute(u,[2 1 3]);
    w_vtk=permute(w,[2 1 3]);
    fwrite(fid, [reshape(u_vtk,1,nr_of_elements);  reshape(v_vtk,1,nr_of_elements); reshape(w_vtk,1,nr_of_elements)],'float','b');
   
    %append scalar data
    fprintf(fid, '\nSCALARS alpha float\n'); %ASCII header
    fprintf(fid, 'LOOKUP_TABLE default\n'); %ASCII header
    alpha_vtk=permute(alpha,[2 1 3]);
    fwrite (fid, reshape(alpha_vtk,1,nr_of_elements),'float','b'); %binary data
    
    fclose(fid);
    
    
end
%% .mat

if handles.saveoptions.mat ==1
    
    filename = [parameter.Name '/' sprintf('%06d',round(ttime*ScaleT,3)) 'Topology'];
    
    save([filename '.mat'],'h');
    
end

%% .fig1

if handles.saveoptions.fig1 ==1
    
    filename = [parameter.Name '/' sprintf('%06d',round(ttime*ScaleT,3)) 'Fig1.fig'];
    orignalAxes = handles.axes2;
    
    Fig1 = figure;
    copyobj(orignalAxes, Fig1);
    pos = get(gcf, 'DefaultAxesPosition'); set(gca, 'Units', 'normalized', 'Position', pos)
    hgsave(Fig1, filename);
    close
end

%% .fig2

if handles.saveoptions.fig2 ==1
    
    if exist([parameter.Name '/' 'Fig2_data.mat']) == 0
        
        %fig2_data(1,1)="time";
        %fig2_data(1,2)="hmin";
        %fig2_data(1,3)="hmax";
        %fig2_data(1,4)="Reynolds / reduced Reynolds";
        fig2_data(2,1)=ttime*ScaleT;
        fig2_data(2,2)=storeData.hmin(end)*ScaleH;
        fig2_data(2,3)=storeData.hmax(end)*ScaleH;
        fig2_data(2,4)=storeData.q(end)*ScaleQ;
        filename = [parameter.Name '/' 'Fig2_data'];
        save([filename '.mat'],'fig2_data');
        
    else
        
        data=load([parameter.Name '/' 'Fig2_data.mat']);
        fig2_data=data.fig2_data;
        fig2_data(size(fig2_data,1)+1,1)=ttime*ScaleT;
        fig2_data(size(fig2_data,1),2)=storeData.hmin(end)*ScaleH;
        fig2_data(size(fig2_data,1),3)=storeData.hmax(end)*ScaleH;
        fig2_data(size(fig2_data,1),4)=storeData.q(end)*ScaleQ;
        filename = [parameter.Name '/' 'Fig2_data'];
        save([filename '.mat'],'fig2_data');
        
    end
    
end

%% Save U

Save_U(parameter,U,ttime,M1,M2, N1, N2, h0, kx, kz, parameter.Shkadov_zeta, parameter.Shkadov_kappa, parameter.DimVar_Lx, parameter.DimVar_Lz, storeData.hmax, ScaleT);

end