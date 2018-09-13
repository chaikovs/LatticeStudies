
%% DIPOLES

data_dip = load('DIPOLE_BI_all.dat');
I = data_dip(:,1);

data_dip9 = load('Cynthia/calib_BI_DIP9.dat');
I9raw = data_dip9(:,1); %current
B9field = data_dip9(:,2); %central field
B9int = data_dip9(:,4); %field integral

data_dip_opera = load('DIPOLE_MagnL_opera.txt');

%%

dipole_meas.dip09full = [data_dip9(:,1) data_dip9(:,4) data_dip9(:,2)];
dipole_meas.dip1to15_100A = [data_dip(4,:)' data_dip(1,:)'];
dipole_meas.dip1to15_200A = [data_dip(5,:)' data_dip(2,:)'];
dipole_meas.dip1to15_275A = [data_dip(6,:)' data_dip(3,:)'];

dipole_meas.lmag_opera = [data_dip_opera(:,1) data_dip_opera(:,2)];

% 
% for idip = 1:15
% 
% dipole_meas.dip1to15{idip} = [data_dip(4:end,idip) data_dip(1:3,idip)];
% 
% end

dipole_meas.dip01 = [data_dip(4:end,1) data_dip(1:3,1)];
dipole_meas.dip02 = [data_dip(4:end,2) data_dip(1:3,2)];
dipole_meas.dip03 = [data_dip(4:end,3) data_dip(1:3,3)];
dipole_meas.dip04 = [data_dip(4:end,4) data_dip(1:3,4)];
dipole_meas.dip05 = [data_dip(4:end,5) data_dip(1:3,5)];
dipole_meas.dip06 = [data_dip(4:end,6) data_dip(1:3,6)];
dipole_meas.dip07 = [data_dip(4:end,7) data_dip(1:3,7)];
dipole_meas.dip08 = [data_dip(4:end,8) data_dip(1:3,8)];
dipole_meas.dip09 = [data_dip(4:end,9) data_dip(1:3,9)];
dipole_meas.dip10 = [data_dip(4:end,10) data_dip(1:3,10)];
dipole_meas.dip11 = [data_dip(4:end,11) data_dip(1:3,11)];
dipole_meas.dip12 = [data_dip(4:end,12) data_dip(1:3,12)];
dipole_meas.dip13 = [data_dip(4:end,13) data_dip(1:3,13)];
dipole_meas.dip14 = [data_dip(4:end,14) data_dip(1:3,14)];
dipole_meas.dip15 = [data_dip(4:end,15) data_dip(1:3,15)];

dipole_meas.fieldIntegral = [data_dip(1:3,1) data_dip(1:3,2) data_dip(1:3,3) data_dip(1:3,4) data_dip(1:3,5)  ...
    data_dip(1:3,6) data_dip(1:3,7) data_dip(1:3,8) data_dip(1:3,9) data_dip(1:3,10) data_dip(1:3,11) data_dip(1:3,12) data_dip(1:3,13) data_dip(1:3,14) data_dip(1:3,15)];

dipole_meas.current = [data_dip(4:end,1) data_dip(4:end,2) data_dip(4:end,3) data_dip(4:end,4) data_dip(4:end,5) ...
    data_dip(4:end,6) data_dip(4:end,7) data_dip(4:end,8) data_dip(4:end,9) data_dip(4:end,10) data_dip(4:end,11) data_dip(4:end,12) data_dip(4:end,13) data_dip(4:end,14) data_dip(4:end,15) ];


save('magn_meas_dipole.mat', 'dipole_meas');

%%


figure
set(gca,'FontSize',18)
plot(dipole_meas.current(1,:), (mean(dipole_meas.fieldIntegral(1,:)) - dipole_meas.fieldIntegral(1,:))./mean(dipole_meas.fieldIntegral(1,:)), 'o', 'MarkerSize',6, 'LineWidth',1.2)
xlabel(' Current [A]')
ylabel('(Mean B - B_i) / Mean B')
title('Spread of the measurements for 15 DIPOLES @100A')


str0= { 'DIP#01' 'DIP#02' 'DIP#03' 'DIP#04' 'DIP#05' 'DIP#06' 'DIP#07' 'DIP#08' 'DIP#09' 'DIP#10' 'DIP#11' 'DIP#12' 'DIP#13' 'DIP#14' 'DIP#15'};
figure
set(gca,'FontSize',18)
plot(dipole_meas.current(1:3,:), (mean(dipole_meas.fieldIntegral(:,:),2) - dipole_meas.fieldIntegral(1:3,:))./mean(dipole_meas.fieldIntegral(:,:),2), 'o-', 'MarkerSize',6, 'LineWidth',1.2)
xlabel(' Current [A]')
ylabel('(Mean B - B_i) / Mean B')
title('Spread of the measurements for 15 DIPOLES')
u = legend(str0);
set(u,'Location','NorthWest')
xlim([30 300])

str= {'DIP#09 (all data)' 'DIP#01' 'DIP#02' 'DIP#03' 'DIP#04' 'DIP#05' 'DIP#06' 'DIP#07' 'DIP#08' 'DIP#09' 'DIP#10' 'DIP#11' 'DIP#12' 'DIP#13' 'DIP#14' 'DIP#15'};
figure
set(gca,'FontSize',20)
plot(dipole_meas.dip09full(:,1), dipole_meas.dip09full(:,2), 'bo-', 'MarkerSize',6, 'LineWidth',1.2)
hold all
plot(dipole_meas.current(1:3,:), dipole_meas.fieldIntegral(:,:), '*', 'MarkerSize',9, 'LineWidth',1.2)
hold off
xlabel(' Current [A]')
ylabel('Inegrated field [T m]')
title('Magnetic calibration for 15 DIPOLES')
% u = legend('DIP-9', 'DIP-91','DIP-7','DIP-11', 'DIP-12');
% set(u,'Location','NorthWest')
u = legend(str);
set(u,'Location','NorthWest')

Lmagn = B9int./B9field;

figure
set(gca,'FontSize',20)
plot(data_dip9(2:end,1), Lmagn(2:end), 'ko-', 'MarkerSize',6, 'LineWidth',1.2)
xlabel(' Current [A]')
ylabel('Magnetic Length [m]')
title('Raw data for the DIPOLE #09')

%% QUADs

close all; clc;
data_quad = load('Qpoles_Int_B2_I.txt');
I = data_quad(:,1);


%%

quad_meas.IntGradient = data_quad(:,2:end);

quad_meas.current = data_quad(:,1);


save('magn_meas_quad.mat', 'quad_meas');


%%

figure(1)
set(gca,'FontSize',18)
plot(I, data_quad(:,2:end), 'o-', 'MarkerSize',8, 'LineWidth',1.2)
xlabel(' Current [A]')
ylabel('Inegrated gradient [T]')
title('Magnetic calibration for 34 QUADS')
% print('quad_all_raw.png','-dpng','-r300')
