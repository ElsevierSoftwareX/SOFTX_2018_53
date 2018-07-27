function [] = Save_U(parameter,U,ttime,M1,M2, N1, N2, h0, kx, kz, zeta, kappa, Lx, Lz,hmax,ScaleT)
% save the solution for a printing time in a mat-file 
Usave=U(end,:);
filename = [parameter.Name '/' sprintf('%06d',round(ttime*ScaleT,3)) '_U_Full2ndOrder.mat'];
save(filename,'Usave', 'M1', 'M2', 'N1', 'N2', 'h0', 'kx', 'kz', 'zeta', 'kappa', 'Lx', 'Lz','ttime','hmax');
end