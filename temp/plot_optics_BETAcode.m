% data1 = load('results_Ring1');

flag_save = 1;

data = load('synoptic_LT_Ring');

data1 = load('results_RingTL1');
data2 = load('results_RingTL2');
data22 = load('results_Ring1');
% data1 = load('results_TL1');
% data2 = load('results_TL2');
%%
figure
plot(data(:,1),data(:,2),'b-', 'Linewidth',2)

%%


figure(22)
set(gca,'FontSize',14)
h1=plot(data22(:,1), 1.5*data22(:,2), 'k-');
axis([0 20 0 5])

%%


figure(2)
set(gca,'FontSize',14)
h1=plot(data1(:,1), 0.5*data1(:,2), 'k-');
hold on
h2=plot(data2(:,1),data2(:,2),'r-', 'Linewidth',1.6);
h3=plot(data2(:,1),data2(:,8),'b-', 'Linewidth',1.6);
h4=plot(data2(:,1),10*data2(:,5),'g-', 'Linewidth',1.6);
hold off
grid on
xlabel('Position [m]')
ylabel('Optical functions [m]');
u = legend([h2 h3 h4],{'\beta_x','\beta_z','10*D'});
%u = legend(' ', '\beta_x','\beta_z','10*D');
set(u,'Location','NorthEast')
axis([0 32 -10 100])
if (flag_save == 1)
    set(gcf, 'color', 'w');
    export_fig(['Optics_TL_Ring'  '.pdf'])
end

%%

%dydx = diff(data2(:,5))./diff(data2(:,1));
dydx = diff([eps; data2(:,5)])./diff([eps; data2(:,1)]);

figure(3)
set(gca,'FontSize',14)
plot(data1(:,1),0.05*data1(:,2), 'k-')
hold on
plot(data2(:,1),data2(:,5),'b-', 'Linewidth',2)
plot(data2(:,1),data2(:,6),'m-', 'Linewidth',2)
%plot(data2(:,1),dydx,'k-')
hold off
xlabel('Position [m]')
ylabel('Beta functions [mm]]');
xlabel('Position [m]')
ylabel('Dispersion functions [mm]]');
u = legend('\eta_x','\eta_x');
set(u,'Location','NorthEast')


%% plot orbit
orb1 = load('orbit1');
orb2 = load('orbit2');
%orb3 = load('orbit3');


figure
set(gca,'FontSize',14)
plot(orb1(:,1),1*orb1(:,2), 'k-')
hold on
plot(orb2(:,1),orb2(:,7),'b-', 'Linewidth',2)
% plot(orb2(:,1),orb2(:,6),'m-', 'Linewidth',2)
% plot(orb2(:,1),orb2(:,7),'m-', 'Linewidth',2)
% plot(orb2(:,1),orb2(:,8),'m-', 'Linewidth',2)
hold off



