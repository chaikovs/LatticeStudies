
close all
clear all


% data_dip9 = load('calib_BI_DIP9.dat');
% I9raw = data_dip9(:,1);
% B9raw = data_dip9(:,4); 
% 
% btemp = (B9raw(30) + B9raw(31))./2;
% itemp = (I9raw(30) + I9raw(31))./2;
% B9 = [B9raw(1:29);btemp; B9raw(32:end)];
% I9 = [I9raw(1:29);itemp; I9raw(32:end)];
% 
% % B9 = B9raw;
% % I9 = I9raw;

load magn_meas_dipole.mat
% I9 = dipole_meas.dip09full(:,1);
% B9 = dipole_meas.dip09full(:,2);

% I9 = [dipole_meas.dip09full(1:30,1); dipole_meas.dip09full(32:end,1)];
% B9 = [dipole_meas.dip09full(1:30,2); dipole_meas.dip09full(32:end,2)];

I9 = [dipole_meas.dip09full(1:29,1); dipole_meas.dip09full(31:end,1)];
B9 = [dipole_meas.dip09full(1:29,2); dipole_meas.dip09full(31:end,2)];

% I9raw = dipole_meas.dip09full(:,1);
% B9raw = dipole_meas.dip09full(:,2);
% btemp = (B9raw(30) + B9raw(31))./2;
% itemp = (I9raw(30) + I9raw(31))./2;
% B9 = [B9raw(1:29);btemp; B9raw(32:end)];
% I9 = [I9raw(1:29);itemp; I9raw(32:end)];


%%

figure(1)
set(gca,'FontSize',16)
plot(I9, B9, 'ko-', 'MarkerSize',6, 'LineWidth',1.2)
% hold on
% plot(I91, B91, 'bd-', 'MarkerSize',6, 'LineWidth',1.2)
% plot(I7, B7, 'ro-', 'MarkerSize',6, 'LineWidth',1.2)
% plot(I11, B12, 'go-', 'MarkerSize',6, 'LineWidth',1.2)
% plot(I12, B12, 'mo-', 'MarkerSize',6, 'LineWidth',1.2)
% hold off
xlabel(' Current [A]')
ylabel('Inegrated field [T m]')
title('Raw data for the DIPOLE #09')
% u = legend('DIP-9', 'DIP-91','DIP-7','DIP-11', 'DIP-12');
% set(u,'Location','NorthWest')
%%

% B = [B9'; B91'; B7'; B11'; B12'];
% Bmean = mean(B);
% Bstd = std(B);
% sizeB = size(B);
% entries = sizeB(1,1);
% 
% table(I9,B9,B91,B7,B11,B12,Bmean',Bstd',Bstd'./Bmean','VariableNames',{'I','B8','B9','B22','Bmean','Bstd','RelError'})

%%

xx9 = linspace(min(I9),max(I9),5*length(I9));
yy9 = spline(I9,B9,xx9);
ii9 = interp1(I9,B9,xx9); 

p9 = polyfit(I9,B9,7);

%x = linspace(min(I8),max(I8),3*length(I8));
f9 = polyval(p9,xx9);

% figure
% set(gca,'FontSize',20)
% plot(I9, B9, 'bo', 'MarkerSize',6)
% hold on
% plot(xx9, yy9, 'r-', 'Linewidth',1.3)
% plot(xx9, ii9, 'k-', 'Linewidth',1.3)
% hold off
% title('DIP 9')
% xlabel(' Current [A]')
% ylabel('Inegrated field [T m]')
% u = legend('Data','Spline', 'Interp1');
% set(u,'Location','NorthWest')


figure(2);
subplot(2,1,1);
set(gca,'FontSize',20)
plot(I9,B9,'-ro',xx9,yy9,'k',xx9,ii9,'b',xx9,f9,'m','MarkerSize',5);
u=legend('Data,','Spline','Interp1','Polynom fit');
set(u,'Location','NorthWest','FontSize',14)
xlabel(' Current [A]')
ylabel('Inegrated field [T m]')
subplot(2,1,2);
set(gca,'FontSize',20)
plot(I9,(spline(xx9,yy9,I9) - B9)./B9,'ro-',...
     I9,(interp1(xx9,ii9,I9) - B9)./B9,'bo-',...
     I9,(polyval(p9,I9) - B9)./B9,'mo-','MarkerSize',5);

u = legend('(Spline - Data)/Data','(Interp1 - Data)/Data','(Polynom fit - Data)/Data');
set(u,'Location','SouthEast','FontSize',14)
xlabel(' Current [A]')
ylabel('Magnetic Field difference')
title('% difference between curve and data');

%% Current range fitting

% B9_reg1 = B9(1:13);
% I9_reg1 = I9(1:13);
% 
% B9_reg2 = B9(13:end);
% I9_reg2 = I9(13:end);
% 
% x1 = linspace(min(I9_reg1),max(I9_reg1),5*length(I9_reg1));
% x2 = linspace(min(I9_reg2),max(I9_reg2),5*length(I9_reg2));
% 
% %%
% 
% p9 = polyfit(I9_reg2,B9_reg2,5)
% 
% f9 = polyval(p9,x2);
% 
% % Linear fit
% pl9 = polyfit(I9_reg1,B9_reg1,1);
% 
% fl9 = polyval(pl9,x1);
% 
% T91 = table(I9_reg1,B9_reg1,polyval(pl9,I9_reg1),B9_reg1-polyval(pl9,I9_reg1),100*(polyval(pl9,I9_reg1) - B9_reg1)./B9_reg1,'VariableNames',{'I','B','Fit','FitError','RelFitError'})
% 
% T92 = table(I9_reg2,B9_reg2,polyval(p9,I9_reg2),B9_reg2-polyval(p9,I9_reg2),100*(polyval(p9,I9_reg2) - B9_reg2)./B9_reg2,'VariableNames',{'I','B','Fit','FitError','RelFitError'})
% 
% 
% figure(3)
% subplot(2,1,1);
% set(gca,'FontSize',16)
% plot(I9, B9, 'ko', 'MarkerSize',6)
% hold on
% plot(x2, f9, 'r-', 'LineWidth',1.3)
% plot(x1, fl9, 'b-', 'LineWidth',1.3)
% hold off
% title('DIP 09 Polynomial fit in 2 regions')
% xlabel(' Current [A]')
% ylabel('Inegrated gradient [T]')
% u = legend('Data','Polynom fit','Linear fit');
% set(u,'Location','NorthWest')
% subplot(2,1,2);
% set(gca,'FontSize',16)
% plot(I9_reg2,(polyval(p9,I9_reg2) - B9_reg2)./B9_reg2,'ro-',...
%      I9_reg1,(polyval(pl9,I9_reg1) - B9_reg1)./B9_reg1,'bo-','MarkerSize',5);
%  xlabel(' Current [A]')
% ylabel('Magnetic Field difference')
% u = legend('(Polynom fit - Data)/Data','(Linear fit - Data)/Data');
% set(u,'Location','SouthEast')


%%

% Linear fit
plinear9 = polyfit(I9,B9,1)

flinear9 = polyval(plinear9,xx9);

figure(4)
set(gca,'FontSize',20)
plot(I9, B9, 'ko', 'MarkerSize',6)
hold on
plot(xx9, flinear9, 'r-', 'LineWidth',1.3)
%plot(x1, fl9, 'b-', 'LineWidth',1.3)
hold off
title('DIPOLE #09 Linear fit')
xlabel(' Current [A]')
ylabel('Inegrated field [T m]')
u = legend('Data','Linear fit');
set(u,'Location','NorthWest','FontSize',14)



%% Current range fitting around working point 

B9_reg1 = B9(1:15);
I9_reg1 = I9(1:15);

B9_reg2 = B9(15:end);
I9_reg2 = I9(15:end);

x1 = linspace(min(I9_reg1),max(I9_reg1),5*length(I9_reg1));
x2 = linspace(min(I9_reg2),max(I9_reg2),5*length(I9_reg2));

%%

p9 = polyfit(I9_reg2,B9_reg2,5)

f9 = polyval(p9,x2);

% Linear fit
pl9 = polyfit(I9_reg1,B9_reg1,1)

fl9 = polyval(pl9,x1);

T91 = table(I9_reg1,B9_reg1,polyval(pl9,I9_reg1),B9_reg1-polyval(pl9,I9_reg1),100*(B9_reg1 - polyval(pl9,I9_reg1))./B9_reg1,'VariableNames',{'I','B','Fit','FitError','RelFitErrorPecent'})

T92 = table(I9_reg2,B9_reg2,polyval(p9,I9_reg2),B9_reg2-polyval(p9,I9_reg2),100*(B9_reg2 - polyval(p9,I9_reg2))./B9_reg2,'VariableNames',{'I','B','Fit','FitError','RelFitErrorPercent'})


figure(5)
subplot(2,1,1);
set(gca,'FontSize',14)
plot(I9, B9, 'ko', 'MarkerSize',5)
hold on
plot(x2, f9, 'r-', 'LineWidth',1.3)
plot(x1, fl9, 'b-', 'LineWidth',1.3)
hold off
title('DIPOLE #09 Polynomial fit')
xlabel(' Current [A]')
ylabel('Inegrated field [T m]')
u = legend('Data','Polynom fit','Linear fit');
set(u,'Location','NorthWest','FontSize',14)
subplot(2,1,2);
set(gca,'FontSize',16)
plot(I9_reg2,(polyval(p9,I9_reg2) - B9_reg2)./B9_reg2,'ro-',...
     I9_reg1,(polyval(pl9,I9_reg1) - B9_reg1)./B9_reg1,'bo-','MarkerSize',5);
 xlabel(' Current [A]')
ylabel('Field difference')
u = legend('(Polynom fit - Data)/Data','(Linear fit - Data)/Data');
set(u,'Location','SouthEast','FontSize',12)
print('dipole_all_fit09.png','-dpng','-r300')

figure(55)
subplot(2,1,1);
set(gca,'FontSize',14)
plot(I9, B9, 'ko', 'MarkerSize',5)
hold on
plot(x2, f9, 'r-', 'LineWidth',1.3)
plot(x1, fl9, 'b-', 'LineWidth',1.3)
hold off
title('DIPOLE #09 Polynomial fit')
xlabel(' Current [A]')
ylabel('Inegrated field [T m]')
u = legend('Data','Polynom fit','Linear fit');
set(u,'Location','NorthWest','FontSize',14)
subplot(2,1,2);
set(gca,'FontSize',16)
plot(I9_reg2,(polyval(p9,I9_reg2) - B9_reg2),'ro-',...
     I9_reg1,(polyval(pl9,I9_reg1) - B9_reg1),'bo-','MarkerSize',5);
 xlabel(' Current [A]')
ylabel('Field difference')
u = legend('Polynom fit - Data','Linear fit - Data');
set(u,'Location','SouthEast','FontSize',12)
print('dipole_all_fit09_simpleDiff.png','-dpng','-r300')


%%

bend2gev(159.2996)
%bend2gev(159.2475)
bend2gev(300)
gev2bend(0.05)

