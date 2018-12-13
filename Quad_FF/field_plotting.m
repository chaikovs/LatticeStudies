
%%

data = xlsread('/Users/ichaikov/Work/ThomX/Bfield_measurements/Quad_Thomx_FF.xlsx')

%%

figure
set(gca,'fontsize',18);
plot(data(:,3), data(:,5),'.-','Markersize',8)
xlabel('z [mm]')
ylabel('By [T] @ R = 18 mm')
title ('x = 18mm , y = 0')
%print('thomx_QUAD_By_FF.png','-dpng','-r300')

figure
plot(data(:,3), data(:,4),'.-')
xlabel('z [mm]')
ylabel('Bx [T]')
title ('x = 18mm, y = 0')


figure
plot(data(:,3), data(:,6),'.-')
xlabel('z [mm]')
ylabel('Bz [T]')
title ('x = 18mm, y = 0')