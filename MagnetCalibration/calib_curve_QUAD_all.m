
close all; clc;
data_quad = load('Qpoles_Int_B2_I.txt');
I = data_quad(:,1);

%%

figure(1)
set(gca,'FontSize',18)
plot(I, data_quad(:,2:end), 'o-', 'MarkerSize',8, 'LineWidth',1.2)
xlabel(' Current [A]')
ylabel('Inegrated gradient [T]')
title('Magnetic calibration for 34 QUADS')
print('quad_all_raw.png','-dpng','-r300')


current_mean = mean(data_quad(:,2:end),2);
current_std = std(data_quad(:,2:end),0,2);

figure(2) 
set(gca,'FontSize',18)
errorbar(I, mean(data_quad(:,2:end),2), std(data_quad(:,2:end),0,2), 'k.-', 'MarkerSize',6, 'LineWidth',1.2)
xlabel(' Current [A]')
ylabel('Inegrated gradient [T]')
title('Magnetic calibration for 34 QUADS (mean value)')


%%

table(I, current_mean, current_std ,current_std./current_mean, 'VariableNames' ,{'I','Bmean','Bstd','RelError'})

%%

xx = linspace(min(I),max(I),5*length(I));
%yy = spline(I,current_mean,xx);
%ii = interp1(I,current_mean,xx); 

% figure
% set(gca,'FontSize',20)
% plot(I, current_mean, 'bo', 'MarkerSize',6)
% hold on
% plot(xx, yy, 'r-', 'Linewidth',1.3)
% plot(xx, ii, 'k-', 'Linewidth',1.3)
% hold off
% title('QUAD mean')
% xlabel(' Current [A]')
% ylabel('Inegrated gradient [T]')
% u = legend('Data','Spline', 'Interp1');
% set(u,'Location','NorthWest')

%% Polynomial fit

p = polyfit(I,current_mean,5)

%x = linspace(min(I8),max(I8),3*length(I8));
f = polyval(p,xx);

%T8 = table(I8,B8,interp1(x,f8,I8),B8-interp1(x,f8,I8),100*(interp1(x,f8,I8) - B8)./B8,'VariableNames',{'I','B','Fit','FitError','RelFitError'})
T = table(I,current_mean,polyval(p,I),current_mean-polyval(p,I),100*(polyval(p,I) - current_mean)./current_mean,'VariableNames',{'I','B','Fit','FitError','RelFitErrorPercent'})

% Linear fit
p1 = polyfit(I,current_mean,1);

f1 = polyval(p1,xx);

T1 = table(I,current_mean,polyval(p1,I),current_mean-polyval(p1,I),100*(polyval(p1,I) - current_mean)./current_mean,'VariableNames',{'I','B','Fit','FitError','RelFitErrorPercent'})


figure(3)
subplot(2,1,1);
set(gca,'FontSize',20)
plot(I, current_mean, 'ko', 'MarkerSize',6)
hold on
plot(xx, f, 'r-', 'LineWidth',1.3)
plot(xx, f1, 'b-', 'LineWidth',1.3)
hold off
title('Magnetic calibration for 34 QUADS (mean value)')
xlabel(' Current [A]')
ylabel('Inegrated gradient [T]')
u = legend('Data','Polynom fit','Linear fit');
set(u,'Location','NorthWest','FontSize',14)
subplot(2,1,2);
set(gca,'FontSize',20)
% plot(I8,(interp1(x,f8,I8) - B8)./B8,'ro-',...
%      I8,(interp1(x,f81,I8) - B8)./B8,'bo-');
plot(I,(polyval(p,I) - current_mean)./current_mean,'ro-',...
     I,(polyval(p1,I) - current_mean)./current_mean,'bo-');
 xlabel(' Current [A]')
ylabel('Magnetic Field difference')
u = legend('(Polynom fit - Data)/Data','(Linear fit - Data)/Data');
set(u,'Location','NorthEast','FontSize',14)
%print('quad_all_fitmean.png','-dpng','-r300')

%% Linear and Polynomial fit

% [fitresultL, gofL] = createFitLin(I(1:end-1), current_mean(1:end-1));
[fitresultL, gofL] = createFitLin(I, current_mean);
[fitresultP, gofP] = createFitPol(I, current_mean);


fdata = feval(fitresultP,xx);

fitparam = coeffvalues(fitresultP)

%% Linear, Polynomial fit and spline

% 
% B_spline = spline(I,current_mean,xx); %pchip
% B_linear = f1;
% 
% figure;
% subplot(2,1,1);
% plot(I,current_mean,'ro',xx, B_spline,'k',xx,B_linear,'b',xx,f,'m','MarkerSize',6);
% u=legend('Data,','Spline','Linear fit','Polynom fit');
% set(u,'Location','NorthWest')
% xlabel(' Current [A]')
% ylabel('Inegrated gradient [T]')
% subplot(2,1,2);
% plot(I,(spline(xx,B_spline,I) - current_mean)./current_mean,'ro-',...
%      I,(polyval(p1,I) - current_mean)./current_mean,'bo-',...
%      I,(polyval(p,I) - current_mean)./current_mean,'mo-');
% u = legend('(Spline - Data)/Data','(Linear fit - Data)/Data','(Polynom fit - Data)/Data');
% set(u,'Location','NorthEast')
% xlabel(' Current [A]')
% ylabel('Magnetic Field difference')
% title('% difference between curve and data');

%%

k2amp('QP2','Monitor',9.2625,[1 2], 0.05)
amp2k('QP2','Monitor',12,[1 2], 0.05)

