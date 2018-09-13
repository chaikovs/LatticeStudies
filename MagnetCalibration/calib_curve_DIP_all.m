

% data_dip = load('DIPOLE_BI_all.dat');
% I = data_dip(:,1);
% 
% data_dip9 = load('calib_BI_DIP9.dat');
% I9raw = data_dip9(:,1);
% B9raw = data_dip9(:,4); 
% 
% btemp = (B9raw(30) + B9raw(31))./2;
% itemp = (I9raw(30) + I9raw(31))./2;
% B9 = [B9raw(1:29);btemp; B9raw(32:end)];
% I9 = [I9raw(1:29);itemp; I9raw(32:end)];

% B9 = B9raw;
% I9 = I9raw;

%%
close all
clear all

load magn_meas_dipole.mat
I9 = [dipole_meas.dip09full(1:29,1); dipole_meas.dip09full(31:end,1)];
B9 = [dipole_meas.dip09full(1:29,2); dipole_meas.dip09full(31:end,2)];

%%

str= {'DIP#09 (all data)' 'DIP#01' 'DIP#02' 'DIP#03' 'DIP#04' 'DIP#05' 'DIP#06' 'DIP#07' 'DIP#08' 'DIP#09' 'DIP#10' 'DIP#11' 'DIP#12' 'DIP#13' 'DIP#14' 'DIP#15'};
figure(1)
set(gca,'FontSize',16)
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
set(u,'Location','SouthEast','FontSize',12)
print('dipole_all_calib.png','-dpng','-r300')

str0= { 'DIP#01' 'DIP#02' 'DIP#03' 'DIP#04' 'DIP#05' 'DIP#06' 'DIP#07' 'DIP#08' 'DIP#09' 'DIP#10' 'DIP#11' 'DIP#12' 'DIP#13' 'DIP#14' 'DIP#15'};
figure(2)
set(gca,'FontSize',16)
plot(dipole_meas.current(1:3,:), (mean(dipole_meas.fieldIntegral(:,:),2) - dipole_meas.fieldIntegral(1:3,:))./mean(dipole_meas.fieldIntegral(:,:),2), 'o-', 'MarkerSize',6, 'LineWidth',1.2)
hold on
plot([159.2475 159.2475 ],[-5e-3 4e-3], 'k-')
hold off
xlabel(' Current [A]')
ylabel('(Mean B - B_i) / Mean B')
title('Spread of the measurements for 15 DIPOLES')
u = legend(str0);
set(u,'Location','NorthWest')
xlim([30 300])
%print('dipole_all_dBB_spread_ALL.png','-dpng','-r300')

%%
% str= {'DIP#09 (all data)' 'DIP#01' 'DIP#02' 'DIP#03' 'DIP#04' 'DIP#05' 'DIP#06' 'DIP#07' 'DIP#08' 'DIP#09' 'DIP#10' 'DIP#11' 'DIP#12' 'DIP#13' 'DIP#14' 'DIP#15'};
% figure(1)
% set(gca,'FontSize',20)
% plot(I9, B9, 'bo-', 'MarkerSize',6, 'LineWidth',1.2)
% hold all
% plot(I, data_dip(:,2:end), '*', 'MarkerSize',9, 'LineWidth',1.2)
% hold off
% xlabel(' Current [A]')
% ylabel('Inegrated field [T m]')
% title('Magnetic calibration for 15 DIPOLES')
% % u = legend('DIP-9', 'DIP-91','DIP-7','DIP-11', 'DIP-12');
% % set(u,'Location','NorthWest')
% u = legend(str);
% set(u,'Location','NorthWest')
% %print('dipole_all_raw.png','-dpng','-r300')

%%
str0= { 'DIP#01' 'DIP#02' 'DIP#03' 'DIP#04' 'DIP#05' 'DIP#06' 'DIP#07' 'DIP#08' 'DIP#09' 'DIP#10' 'DIP#11' 'DIP#12' 'DIP#13' 'DIP#14' 'DIP#15'};
figure(3)
set(gca,'FontSize',16)
%plot(I, data_dip(:,2:end), 'o-', 'MarkerSize',9, 'LineWidth',1.2)
plot(dipole_meas.current(1:3,:), dipole_meas.fieldIntegral(:,:), '*-', 'MarkerSize',9, 'LineWidth',1.2)
xlabel(' Current [A]')
ylabel('Inegrated field [T m]')
title('Magnetic calibration for 15 DIPOLES')
u = legend(str0);
set(u,'Location','SouthEast','FontSize',12)
print('dipole_all3P_raw.png','-dpng','-r300')

%%

% figure(22)
% set(gca,'FontSize',18)
% plot(I(1), mean(data_dip(1,2:end),2)-data_dip(1,2:end), 'o-', 'MarkerSize',6, 'LineWidth',1.2)
% hold on
% plot(I(2), mean(data_dip(2,2:end),2)-data_dip(2,2:end), 'o-', 'MarkerSize',6, 'LineWidth',1.2)
% plot(I(3), mean(data_dip(3,2:end),2)-data_dip(3,2:end), 'o-', 'MarkerSize',6, 'LineWidth',1.2)
% hold off
% xlabel(' Current [A]')
% ylabel('Mean field B - B_i [T m]')
% title('Spread of the measurements for 15 DIPOLES')
% u = legend(str0);
% set(u,'Location','NorthWest')
% xlim([30 300])
% %print('dipole_all_measSpread.png','-dpng','-r300')
% 
% figure(222)
% set(gca,'FontSize',18)
% plot(I(1), (mean(data_dip(1,2:end),2)-data_dip(1,2:end).^2./mean(data_dip(1,2:end),2)), 'o-', 'MarkerSize',6, 'LineWidth',1.2)
% hold on
% plot(I(2), (mean(data_dip(2,2:end),2)-data_dip(2,2:end).^2./mean(data_dip(2,2:end),2)), 'o-', 'MarkerSize',6, 'LineWidth',1.2)
% plot(I(3), (mean(data_dip(3,2:end),2)-data_dip(3,2:end).^2./mean(data_dip(3,2:end),2)), 'o-', 'MarkerSize',6, 'LineWidth',1.2)
% hold off
% xlabel(' Current [A]')
% ylabel('Mean field B - B_i / Mean field')
% title('Spread of the measurements for 15 DIPOLES')
% u = legend(str0);
% set(u,'Location','NorthWest')
% xlim([30 300])
% %print('dipole_all_measSpreadRel.png','-dpng','-r300')

% figure
% set(gca,'FontSize',18)
% plot(I, mean(data_dip(:,2:end),2)-data_dip(:,2:end), 'o-', 'MarkerSize',6, 'LineWidth',1.2)
% xlabel(' Current [A]')
% ylabel('Mean B - B_i [T m]')
% title('Spread of the measurements for 15 DIPOLES')
% u = legend(str0);
% set(u,'Location','NorthWest')
% xlim([30 300])
% print('dipole_all_measSpreadJOIN.png','-dpng','-r300')
% 
% figure
% set(gca,'FontSize',18)
% plot(I, (mean(data_dip(:,2:end),2)-data_dip(:,2:end))./mean(data_dip(:,2:end),2), 'o-', 'MarkerSize',6, 'LineWidth',1.2)
% xlabel(' Current [A]')
% ylabel('(Mean B - B_i) / Mean B')
% title('Spread of the measurements for 15 DIPOLES')
% u = legend(str0);
% set(u,'Location','NorthWest')
% xlim([30 300])
% print('dipole_all_measSpreadRelJOIN.png','-dpng','-r300')


%%

% current_mean = mean(data_dip(:,2:end),2);
% current_std = std(data_dip(:,2:end),0,2);
% 
% figure(3)
% set(gca,'FontSize',18)
% errorbar(I, mean(data_dip(:,2:end),2), std(data_dip(:,2:end),0,2), '.-', 'MarkerSize',6, 'LineWidth',1.2)
% xlabel(' Current [A]')
% ylabel('Inegrated field [T m]')
% title('Magnetic calibration for 15 DIPOLES')


%%

% table(I, current_mean, current_std ,current_std./current_mean, 'VariableNames' ,{'I','Bmean','Bstd','RelError'})

%% Preparation for the sorting

% dBBmeas_temp = zeros(4,15);
% dBBmeas_temp(1,:) = 1:15;
% dBBmeas_temp(2:end,:) = (mean(data_dip(:,2:end),2)-data_dip(:,2:end))./mean(data_dip(:,2:end),2);
% dBBmeas = [dBBmeas_temp(:,2) dBBmeas_temp(:,4:end)];
% dBBmeas_mean = mean(dBBmeas(2:end,:),2);
% dBBmeas_std = std(dBBmeas(2:end,:),0,2);
% 
% str_sort = { 'DIP#02' 'DIP#04' 'DIP#05' 'DIP#06' 'DIP#07' 'DIP#08' 'DIP#09' 'DIP#10' 'DIP#11' 'DIP#12' 'DIP#13' 'DIP#14' 'DIP#15'};
% 
% figure
% set(gca,'FontSize',18)
% plot(I, dBBmeas(2:end,:), 'o-', 'MarkerSize',6, 'LineWidth',1.2)
% xlabel(' Current [A]')
% ylabel('(Mean B - B_i) / Mean B')
% title('Spread of the measurements for 15 DIPOLES')
% u = legend(str_sort);
% set(u,'Location','NorthWest')
% xlim([30 300])
% print('dipole_all_measSpreadRelJOIN_sort.png','-dpng','-r300')

