
clear all; close all; clc
%%

load dipole_09_160A_ALBA.dat
load ThomX_dipole_trajectoryALL160A.mat

% load ThomX_dipole_traj160A.mat
% load ThomX_dipole_trajGEOM160A.mat

%%

s_meas = dipole_09_160A_ALBA(:,1);
x_meas = dipole_09_160A_ALBA(:,3);
z_meas = dipole_09_160A_ALBA(:,2);
B1_meas = dipole_09_160A_ALBA(:,6);

%%

figure
set(gcf,'color','w')
set(gca,'fontsize',18)
plot(1e3.*s_meas, B1_meas, 'k.-','DisplayName', 'Meas ALBA 160 A' )
hold on
plot(Sr, Bztrack, 'r.-','DisplayName', 'Tracking meas 160 A' )
plot(1e3.*S0, Bz0, 'b.-','DisplayName', 'Ideal path meas 160 A' )
%plot(1e3.*s_meas200A, B1_meas200A, 'b.-','DisplayName', 'Meas ALBA 200 A Christelle' )
%plot(1e3.*s_simu, B1_simu, 'r.-','DisplayName', 'Simu OPERA 200 A' )
%plot(1e3.*xq, vq, 'g*-')
hold off
xlabel('s-traj [mm]');
ylabel('B1@Rref [T]');
u = legend('show','Location','SouthWest');
set(u,'FontSize',14)
print('DIPOEL_B1_ALBA_geom_tracking_160A','-dpng','-r300')

%%

sinterp = linspace(0,290,400);
B_meas_interp = interp1(1e3.*s_meas, B1_meas,sinterp);
Br_interp = interp1(Sr, Bztrack,sinterp);
B0_interp = interp1(1e3.*S0, Bz0,sinterp);


% figure
% plot(sinterp,B_meas_interp,'bo')
% hold on
% plot(1e3*s_meas, B1_meas, 'r-' )
% hold off

dB1 = B_meas_interp - Br_interp;
dB2 = B_meas_interp - B0_interp;
figure(3)
set(gcf,'color','w')
set(gca,'fontsize',18)
plot(sinterp,dB1, 'r.-','DisplayName', 'Diff ALBA - Tracking' )
hold on
plot(sinterp,dB2, 'k.-','DisplayName', 'Diff ALBA - Geom Traj' )
hold off
xlim([0 max(sinterp)])
xlabel('s-traj (mm)')
ylabel('dB (T)')
grid on
u = legend('show','Location','SouthWest');
set(u,'FontSize',14)
print('DIPOEL_160A_Bs_diff','-dpng','-r300')
