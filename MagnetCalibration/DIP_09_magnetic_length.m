
%%
close all
% data_dip9 = load('calib_BI_DIP9.dat');
% I9 = data_dip9(2:end,1);
% B9field = data_dip9(2:end,2);
% B9int = data_dip9(2:end,4); 
% 
% Lmagn = B9int./B9field;

load magn_meas_dipole.mat


I9 = [dipole_meas.dip09full(1:29,1); dipole_meas.dip09full(31:end,1)];
B9int = [dipole_meas.dip09full(1:29,2); dipole_meas.dip09full(31:end,2)];
B9field = [dipole_meas.dip09full(1:29,3); dipole_meas.dip09full(31:end,3)];

Lmagn = B9int./B9field;

%%

figure(1)
set(gca,'FontSize',16)
plot(I9, B9field, 'ko-', 'MarkerSize',6, 'LineWidth',1.2)
xlabel(' Current [A]')
ylabel('Magnetic field [T]')
title('Raw data for the DIPOLE #09')



figure(2)
set(gca,'FontSize',16)
plot(I9, B9int, 'ko-', 'MarkerSize',6, 'LineWidth',1.2)
xlabel(' Current [A]')
ylabel('Inegrated field [T m]')
title('Raw data for the DIPOLE #09')



figure(3)
set(gca,'FontSize',16)
plot(I9(2:end), Lmagn(2:end), 'ko-', 'MarkerSize',6, 'LineWidth',1.2)
xlabel(' Current [A]')
ylabel('Magnetic Length [m]')
title('Raw data for the DIPOLE #09')


figure(4)
set(gca,'FontSize',20)
plot(I9(2:end), Lmagn(2:end), 'ko-', 'MarkerSize',6, 'LineWidth',1.2,'DisplayName', 'DIP#09 data')
hold on
plot(dipole_meas.lmag_opera(:,2), 1e-3*dipole_meas.lmag_opera(:,1), 'ro-', 'MarkerSize',6, 'LineWidth',1.2,'DisplayName', 'Simulations')
hold off
legend('show','Location','NorthEast');
xlabel(' Current [A]')
ylabel('Magnetic Length [m]')
title('Raw data for the DIPOLE #09 vs Simulations')
print('dipole09_Lmag_data_simu','-dpng','-r300')
%Simu @159 A Lmag = 297.1243 mm, Meas at 159.2996 A Lmag = 296.2326 mm
%=>>0.3% diff




%%

%xx9 = linspace(min(I9),max(I9),5*length(I9));

BF_reg1 = B9field(1:15);
BI_reg1 = B9int(1:15);
I9_reg1 = I9(1:15);

BF_reg2 = B9field(15:end);
BI_reg2 = B9int(15:end);
I9_reg2 = I9(15:end);

x1 = linspace(min(I9_reg1),max(I9_reg1),5*length(I9_reg1));
x2 = linspace(min(I9_reg2),max(I9_reg2),5*length(I9_reg2));

%%
% bend2gev(159.2475)

Inom = gev2bend(0.05); %  159.29959 A

pF9_lin = polyfit(I9_reg1,BF_reg1,2);
pI9_lin = polyfit(I9_reg1,BI_reg1,2);

pF9 = polyfit(I9_reg2,BF_reg2,5);
pI9 = polyfit(I9_reg2,BI_reg2,5);

figure(5)
set(gca,'FontSize',20)
plot(I9_reg1, BI_reg1,'bo', I9_reg2, BI_reg2,'ro' , 'MarkerSize',6, 'LineWidth',1.2)
hold on
plot(x1, polyval(pI9_lin,x1),'k-', x2, polyval(pI9,x2),'k-', 'MarkerSize',6, 'LineWidth',1.2)
hold off


figure(6)
set(gca,'FontSize',20)
plot(I9_reg1, BF_reg1,'bo', I9_reg2,BF_reg2,'ro' , 'MarkerSize',6, 'LineWidth',1.2)
hold on
plot(x1, polyval(pF9_lin,x1),'k-', x2, polyval(pF9,x2),'k-', 'MarkerSize',6, 'LineWidth',1.2)
hold off

%%
%magn_length_nominal = polyval(pI9,159.2475)./polyval(pF9,159.2475);

magn_length_nominal = polyval(pI9,Inom)./polyval(pF9,Inom)

magn_length_160A = polyval(pI9,160)./polyval(pF9,160)

