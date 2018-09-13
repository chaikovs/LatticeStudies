clear all; close all;

%all_errors = load('data_ALLerrorsDA_align_err_30um_02mrad_field_err_0001');
all_errors = load('data_ALLerrorsMultipBPM_align_err_100um_05mrad_field_err_0005_dpp_0_GOOD');

%%
figure(11);
plot(1e3*all_errors.E.xorbitmax,1e3*all_errors.E.yorbitmax,'+r','MarkerSize',8)
set(gcf,'color','w')
set(gca,'fontsize',16');
xlim([0 25])
ylim([0 15])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x_{max} [mm]');
ylabel('y_{max} [mm]');
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLmaxErrors_orbitmax.png','-dpng','-r300')


figure(12);
plot(1e3*all_errors.E.xorbitrms,1e3*all_errors.E.yorbitrms,'xr','MarkerSize',8)
set(gcf,'color','w')
set(gca,'fontsize',16');
xlim([0 15])
ylim([0 10])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x_{rms} [mm]');
ylabel('y_{rms} [mm]');
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLmaxErrors_orbitrms.png','-dpng','-r300')


dbetax=squeeze(all_errors.E.dbeta(:,1,:));sdbetax=std(dbetax');maxdbetax=max(dbetax');dbetaxmax=max(abs(dbetax));
dbetaz=squeeze(all_errors.E.dbeta(:,2,:));sdbetaz=std(dbetaz');maxdbetaz=max(dbetaz');dbetazmax=max(abs(dbetaz));

figure(13);
plot(all_errors.E.SPos,sdbetax*100,'-r', 'Linewidth', 1.6);hold on
plot(all_errors.E.SPos,sdbetaz*100,'-b', 'Linewidth', 1.6);hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
legend('Horizontal','Vertical')
xlabel('s-position [m]');               
ylabel('[\Delta\beta/\beta]_{rms} [%]')
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLmaxErrors_betabeatrms_vs_s.png','-dpng','-r300')

figure(14);
plot(all_errors.E.SPos,maxdbetax*100,'-r', 'Linewidth', 1.6);hold on
plot(all_errors.E.SPos,maxdbetaz*100,'-b', 'Linewidth', 1.6);hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
legend('Horizontal','Vertical')
xlabel('s-position [m]');               
ylabel('[\Delta\beta/\beta]_{max} [%]')
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLmaxErrors_betabeatmax_vs_s.png','-dpng','-r300')

figure(15);
plot(dbetaxmax,dbetazmax,'ro');
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('[\Delta\beta_x/\beta_x]_{max}');
ylabel('[\Delta\beta_z/\beta_z]_{max}');
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLmaxErrors_betabeatmax.png','-dpng','-r300')

figure(16);
plot(all_errors.E.tunes(:,1),all_errors.E.tunes(:,2),'ob','MarkerSize',8);
hold on
plot(all_errors.E.tunes0(:,1),all_errors.E.tunes0(:,2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',9)
hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
u=legend('Random seeds','Initial tunes (3.17/1.64)')
set(u, 'Location','NorthEast')
xlabel('\nu_x');
ylabel('\nu_z');
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLmaxErrors_tunes.png','-dpng','-r300')


%% DA

figure(17);
plot(all_errors.E.xda0*1e3,all_errors.E.zda0*1e3,'b-','LineWidth',2);hold on,
plot(all_errors.E.xda*1e3,all_errors.E.zda*1e3,'color',[0,0,0]+0.5,'LineWidth',1);
plot(all_errors.E.xda0*1e3,all_errors.E.zda0*1e3,'b-','LineWidth',2);hold off
legend('Initial: no errors','Trials with errors')
xlim([-50 50])
ylim([0 30])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x [mm]');                 % Add labels
ylabel('z [mm]');
title('DA (\deltap = 0)')
%print('thomx_ALLmaxErrors_DA_init_seeds.png','-dpng','-r300')

% Stat over DA
xxdai=[-60:0.2:60]*1e-3;
zzdai0=interp1(all_errors.E.xda0,all_errors.E.zda0,xxdai','pchip',0);
for kerr=1:486%all_errors.E.nerr
    zzdai(:,kerr)=interp1(all_errors.E.xda(:,kerr),all_errors.E.zda(:,kerr),xxdai','pchip',0);
end

mzzdai=mean(zzdai,2);
szzdai=std(zzdai')';
DA_surf=sum(zzdai)*0.1e-3; % surface in m^2
DA_surf0=sum(zzdai0)*0.1e-3; % surface in m^2

figure(18)
plot(xxdai*1e3,zzdai0*1e3,'-b','LineWidth',2);hold on
plot(xxdai*1e3,mzzdai*1e3,'-r','LineWidth',2);
plot(xxdai*1e3,(mzzdai-szzdai)*1e3,'--m','LineWidth',2);
plot(xxdai*1e3,(mzzdai+szzdai)*1e3,'--m','LineWidth',2);hold off
legend('Initial: no errors','Mean of 500 seeds','Mean - \sigma','Mean + \sigma')
xlim([-60 60])
ylim([0 25])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x [mm]');                 % Add labels
ylabel('z [mm]');
title('DA (\deltap = 0)')
%print('thomx_ALLmaxErrors_DA_init_,meanseeds1.png','-dpng','-r300')

figure(181)
plot(all_errors.E.xda0*1e3,all_errors.E.zda0*1e3,'b-','LineWidth',2);hold on
plot(xxdai*1e3,mzzdai*1e3,'-r','LineWidth',2);
plot(all_errors.E.xda*1e3,all_errors.E.zda*1e3,'color',[0,0,0]+0.5,'LineWidth',1);
plot(xxdai*1e3,mzzdai*1e3,'-r','LineWidth',2);
plot(all_errors.E.xda0*1e3,all_errors.E.zda0*1e3,'b-','LineWidth',2)
% plot(xxda0*1e3,zzda0*1e3,'m-','LineWidth',2)
%plot(xxdai*1e3,(mzzdai-szzdai)*1e3,'--m','LineWidth',2);
%plot(xxdai*1e3,(mzzdai+szzdai)*1e3,'--m','LineWidth',2);hold off
legend('Initial: no errors','Mean of 500 seeds')
xlim([-60 60])
ylim([0 30])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x [mm]');                 % Add labels
ylabel('z [mm]');
title('DA (\deltap = 0)')
%print('thomx_ALLmaxErrors_DA_init_meanseeds2.png','-dpng','-r300')

figure(19);
plot(DA_surf0*1e6, 0 ,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10);hold on
h=histogram(DA_surf*1e6); 
h.FaceColor = [0 0.5 0.5];
hold off    
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('DA surface [mm^2]');
ylabel('Entries');
u=legend('Initial: no errors','Trials with errors')
set(u,'Location','NorthWest','Orientation','vertical')
%print('thomx_ALLmaxErrors_DA_surf.png','-dpng','-r300')

%%

% figure('units','normalized','position',[0.3 0.3 0.45 0.35])
% atplot(rerr,'comment',[],@plClosedOrbit)



figure('units','normalized','position',[0.1 0.4 0.45 0.35])
atplot(all_errors.E.ringerrors(:,1),'comment',[],@pltmisalignments);
%print('thomx_ALLmaxErrors_misalignments.png','-dpng','-r300')





