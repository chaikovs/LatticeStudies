
clear all; close all; clc
%%

load dipole_09_160A_ALBA.dat
load ThomX_dipole_trajectoryALL160A.mat

%%

s_meas = dipole_09_160A_ALBA(:,1);
x_meas = dipole_09_160A_ALBA(:,3);
z_meas = dipole_09_160A_ALBA(:,2);
B1_meas = dipole_09_160A_ALBA(:,6);

%%

figure
set(gcf,'color','w')
set(gca,'fontsize',18)
plot(1e3*z_meas, 1e3*x_meas, '.-','DisplayName', 'Meas ALBA 160 A' )
hold on
plot(Yr, Xr, 'r.-','DisplayName', 'Tracking in Bfield 160 A' )
%plot(1e3*S0, X0, 'k.-','DisplayName', 'Geometric traj 160 A' )
plot(Y0, X0, 'k.-','DisplayName', 'Geometric trajectory 160 A' )
hold off
xlabel('s [mm]');
ylabel('x [mm]');
u = legend('show','Location','SouthWest');
set(u,'FontSize',14)
print('DIPOEL_160A_traj_comp','-dpng','-r300')

yinterp = linspace(0,300,400);
x_meas_interp = interp1(1e3*z_meas, 1e3*x_meas,yinterp);
xr_interp = interp1(Yr, Xr,yinterp);
x0_interp = interp1(Y0, X0,yinterp);


% figure
% plot(sinterp,x_meas_interp,'bo')
% hold on
% plot(1e3*s_meas, 1e3*x_meas, 'r-' )
% hold off

dx1 = x_meas_interp - xr_interp;
dx2 = x_meas_interp - x0_interp;
figure(3)
set(gcf,'color','w')
set(gca,'fontsize',18)
plot(yinterp,dx1, 'r.-','DisplayName', 'Diff ALBA - Tracking' )
hold on
plot(yinterp,dx2, 'k.-','DisplayName', 'Diff ALBA - Geom Traj' )
hold off
xlim([0 max(yinterp)])
xlabel('s (mm)')
ylabel('dx (mm)')
grid on
u = legend('show','Location','SouthWest');
set(u,'FontSize',14)
print('DIPOEL_160A_traj_diff','-dpng','-r300')

%%
90+atan((S0(end)-S0(end-1)).*1e3./(X0(end)-X0(end-1)))*180/pi

90+atan((Y0(end)-Y0(end-1))./(X0(end)-X0(end-1)))*180/pi

%%
% 
% load ThomX_dipole_trajectoryALL_5mm.mat
% figure
% plot(S, Xr, 'r-','DisplayName', 'Tracking in Bfield 160 A' )
% hold on
% plot(1e3*S0, X0, 'k-','DisplayName', 'Geometric traj 160 A' )
% load ThomX_dipole_trajectoryALL_min5mm.mat
% plot(S, Xr, 'r-','DisplayName', 'Tracking in Bfield 160 A' )
% plot(1e3*S0, X0, 'k-','DisplayName', 'Geometric traj 160 A' )
% load ThomX_dipole_trajectoryALL_min10mm.mat
% plot(S, Xr, 'r-','DisplayName', 'Tracking in Bfield 160 A' )
% plot(1e3*S0, X0, 'k-','DisplayName', 'Geometric traj 160 A' )
% load ThomX_dipole_trajectoryALL_10mm.mat
% plot(S, Xr, 'r-','DisplayName', 'Tracking in Bfield 160 A' )
% plot(1e3*S0, X0, 'k-','DisplayName', 'Geometric traj 160 A' )
% load ThomX_dipole_trajectoryALL_min20mm.mat
% plot(S, Xr, 'r-','DisplayName', 'Tracking in Bfield 160 A' )
% plot(1e3*S0, X0, 'k-','DisplayName', 'Geometric traj 160 A' )
% load ThomX_dipole_trajectoryALL_20mm.mat
% plot(S, Xr, 'r-','DisplayName', 'Tracking in Bfield 160 A' )
% plot(1e3*S0, X0, 'k-','DisplayName', 'Geometric traj 160 A' )
% load ThomX_dipole_trajectoryALL_min15mm.mat
% plot(S, Xr, 'r-','DisplayName', 'Tracking in Bfield 160 A' )
% plot(1e3*S0, X0, 'k-','DisplayName', 'Geometric traj 160 A' )
% load ThomX_dipole_trajectoryALL_15mm.mat
% plot(S, Xr, 'r-','DisplayName', 'Tracking in Bfield 160 A' )
% plot(1e3*S0, X0, 'k-','DisplayName', 'Geometric traj 160 A' )

%%
[Xb,Sb,Bz] = getMEASUREDfield2D('ThomX-dipole/field/20160208a_THOMX#009_fieldmap_160A.xlsx');

load ThomX_dipole_trajectoryALL160A.mat
figure
set(gcf,'color','w')
set(gca,'fontsize',18)
plot(Y0, X0,'HandleVisibility','off' )
hold on
imagesc(Sb(:,1),Xb(1,:),Bz');
shading interp;
xlim([-300 300])
ylim([-100 30])
plot(Y0, X0, 'k-','DisplayName', 'Ideal path 160 A Ref' )
plot(Yr, Xr, 'r-','DisplayName', 'Tracking in Bfield 160 A Ref' )
load ThomX_dipole_trajectoryALL160A_p20mm.mat
plot(Yr, Xr, 'r--','DisplayName', 'Tracking in Bfield 160 A +20 mm' )
plot(Y0, X0, 'k--','DisplayName', 'Ideal path 160 A +20 mm' )
load ThomX_dipole_trajectoryALL160A_m20mm.mat
plot(Yr, Xr, 'r-.','DisplayName', 'Tracking in Bfield 160 A -20 mm' )
plot(Y0, X0, 'k-.','DisplayName', 'Ideal path 160 A -20 mm' )
u = legend('show','Location','SouthWest');
set(u,'FontSize',14)
xlabel('s (mm)')
ylabel('x (mm)')
print('DIPOEL_160A_track_3traj','-dpng','-r300')


%%
figure
plot(Yr,'.-')
hold on
plot(Sr,'.-')




