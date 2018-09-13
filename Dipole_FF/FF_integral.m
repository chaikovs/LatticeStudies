
clear all; close all; clc
%%

load dipole_09_160A_ALBA.dat
load ThomX_dipole_trajectoryALL160A.mat
% load ThomX_dipole_traj160A.mat
% load ThomX_dipole_trajGEOM160A.mat

% data_meas200A = load('profilBz_mesure200A_225deg_rho373_corrlin.txt');
% data_simu = load('profilBz_operaV33XC10_225deg_rho373_corrlin.txt');

%%

% s_meas = dipole_09_200A_ALBA(:,1);
% x_meas = dipole_09_200A_ALBA(:,3);
% z_meas = dipole_09_200A_ALBA(:,2);
%B1_meas = dipole_09_200A_ALBA(:,6);

s_meas = dipole_09_160A_ALBA(:,1);
x_meas = dipole_09_160A_ALBA(:,3);
z_meas = dipole_09_160A_ALBA(:,2);
B1_meas = dipole_09_160A_ALBA(:,6);%.*0.996;

% s_meas200A =  1e-3.*data_meas200A(:,1);
% B1_meas200A = data_meas200A(:,2);
% 
% s_simu = 1e-3.*data_simu(:,1);
% B1_simu = data_simu(:,2);

gap = 0.042;

%%

figure
plot(1e3.*s_meas, B1_meas, 'k.-','DisplayName', 'Meas ALBA 160 A' )
hold on
plot(Sr, Bztrack, 'r.-','DisplayName', 'Tracking meas 160 A' )
plot(1e3.*S0, Bz0, 'b.-','DisplayName', 'Ideal path meas 160 A' )
%plot(1e3.*s_meas200A, B1_meas200A, 'b.-','DisplayName', 'Meas ALBA 200 A Christelle' )
%plot(1e3.*s_simu, B1_simu, 'r.-','DisplayName', 'Simu OPERA 200 A' )
%plot(1e3.*xq, vq, 'g*-')
hold off
xlabel('s [mm]');
ylabel('B1@Rref [T]');
u = legend('show','Location','SouthWest');
set(u,'FontSize',14)
%print('DIPOEL_B1_simu_meas','-dpng','-r300')

% figure
% plot( B1_meas, 'k.-','DisplayName', 'Meas ALBA 200 A' )
% hold on
% plot(B1_meas200A, 'b.-','DisplayName', 'Meas ALBA 200 A Christelle' )
% plot(B1_simu, 'r.-','DisplayName', 'Simu OPERA 200 A' )
% %plot(1e3.*xq, vq, 'g*-')
% hold off
% xlabel('Point number');
% ylabel('B1@Rref [T]');
% u = legend('show','Location','NorthEast');
% set(u,'FontSize',14)
% print('DIPOEL_B1_simu_meas2','-dpng','-r300')


%% interp Inom


% xq = -0.050:0.005:0.050;
% vq = interp1(s_meas200A(43:70),B1_meas200A(43:70),xq);
% B0_meas200A = interp1(s_meas200A(43:70),B1_meas200A(43:70),0);
% B0_simu = interp1(s_simu(43:70),B1_simu(43:70),0);
%     
% 
% 
% B1int_meas = trapz(s_meas,B1_meas);
% B1int_meas200A = trapz(s_meas200A,B1_meas200A);
% B1int_simu = trapz(s_simu,B1_simu);
% 
% L_meas = B1int_meas/B1_meas(find(s_meas==0))
% L_meas200A = B1int_meas200A/B0_meas200A
% L_simu = B1int_simu/B0_simu

%%

K_meas = 1./(gap.*B1_meas(find(s_meas==0)).^2) * trapz(s_meas,(B1_meas.*(B1_meas(find(s_meas==0)) - B1_meas)))./2
%K_meas200A = 1./(gap.*B0_meas200A.^2) * trapz(s_meas200A,(B1_meas200A.*(B0_meas200A - B1_meas200A)))./2
%K_simu = 1./(gap.*B0_simu.^2) * trapz(s_simu,(B1_simu.*(B0_simu - B1_simu)))./2

%%

figure
plot(s_meas, x_meas, '.-')

figure
plot(s_meas, z_meas, '.-')


figure
plot(x_meas,s_meas, '.-')
