

flag_save = 0;

 data1 = load('optics_13032017');
 data2 = load('synoptic');

%%
figure
plot(data2(:,1),data2(:,2),'b-', 'Linewidth',2)

%%



figure(2)
set(gca,'FontSize',14)
h1=plot(data2(:,1), 0.5*data2(:,2), 'k-');
hold on
h2=plot(data1(:,1),data1(:,2),'r-', 'Linewidth',1.6);
h3=plot(data1(:,1),data1(:,8),'b-', 'Linewidth',1.6);
h4=plot(data1(:,1),10*data1(:,5),'g-', 'Linewidth',1.6);
hold off
grid on
xlabel('Position [m]')
ylabel('Optical functions [m]');
u = legend([h2 h3 h4],{'\beta_x','\beta_z','10*D'});
set(u,'Location','NorthEast')
if (flag_save == 1)
    set(gcf, 'color', 'w');
    export_fig(['Optics_Ring'  '.pdf'])
end

%%







