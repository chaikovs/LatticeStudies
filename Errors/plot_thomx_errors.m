clear all; close all;

% % error max
% DIP_align = load('data_dipole_align_err_100um_05mrad');
% DIP_field = load('data_dipole_field_err_0005');
% QUAD_align = load('data_quad_align_err_100um_05mrad');
% QUAD_field = load('data_quad_field_err_0005');
% SEXT_align = load('data_sext_align_err_100um_05mrad');
% SEXT_field = load('data_sext_field_err_0005');

% error min
DIP_align = load('data_dipole_align_err_30um_02mrad');
DIP_field = load('data_dipole_field_err_0001');
QUAD_align = load('data_quad_align_err_30um_02mrad');
QUAD_field = load('data_quad_field_err_0001');
SEXT_align = load('data_sext_align_err_30um_02mrad');
SEXT_field = load('data_sext_field_err_0001');

%%


figure(11);
%plot(xorbitmax,yorbitmax,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10);hold on
plot(1e3*DIP_align.E.xorbitmax,1e3*DIP_align.E.yorbitmax,'+b','MarkerSize',8,'DisplayName', 'Dipole alignment errors')
hold on
plot(1e3*DIP_field.E.xorbitmax,1e3*DIP_field.E.yorbitmax,'+r','MarkerSize',8,'DisplayName','Dipole field errors')
plot(1e3*QUAD_align.E.xorbitmax,1e3*QUAD_align.E.yorbitmax,'+m','MarkerSize',8,'DisplayName', 'Quadrupole alignment errors')
plot(1e3*QUAD_field.E.xorbitmax,1e3*QUAD_field.E.yorbitmax,'+g','MarkerSize',8,'DisplayName','Quadrupole field errors')
plot(1e3*SEXT_align.E.xorbitmax,1e3*SEXT_align.E.yorbitmax,'+c','MarkerSize',8,'DisplayName', 'Sextupole alignment errors')
%plot(1e3*SEXT_field.E.xorbitmax,1e3*SEXT_field.E.yorbitmax,'+b')
hold off
set(gcf,'color','w')
set(gca,'fontsize',16');
% xlim([-0.5 20])
% ylim([-0.5 10])
xlim([-0.1 6])
ylim([-0.1 3])
grid on
xlabel('x_{max} [mm]');
ylabel('y_{max} [mm]');
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
print('thomx_errors_orbitmax.png','-dpng','-r300')


figure(12);
%plot(xorbitmax,yorbitmax,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10);hold on
plot(1e3*DIP_align.E.xorbitrms,1e3*DIP_align.E.yorbitrms,'xb','MarkerSize',9,'DisplayName', 'Dipole alignment errors')
hold on
plot(1e3*DIP_field.E.xorbitrms,1e3*DIP_field.E.yorbitrms,'xr','MarkerSize',9,'DisplayName','Dipole field errors')
plot(1e3*QUAD_align.E.xorbitrms,1e3*QUAD_align.E.yorbitrms,'xm','MarkerSize',9,'DisplayName', 'Quadrupole alignment errors')
plot(1e3*QUAD_field.E.xorbitrms,1e3*QUAD_field.E.yorbitrms,'xg','MarkerSize',9,'DisplayName','Quadrupole field errors')
plot(1e3*SEXT_align.E.xorbitrms,1e3*SEXT_align.E.yorbitrms,'xc','MarkerSize',9,'DisplayName', 'Sextupole alignment errors')
%plot(1e3*SEXT_field.E.xorbitrms,1e3*SEXT_field.E.yorbitrms,'xg','MarkerSize',9)
hold off
set(gcf,'color','w')
set(gca,'fontsize',16');
% xlim([-0.5 10])
% ylim([-0.5 6])
xlim([-0.05 3])
ylim([-0.05 2])
grid on
xlabel('x_{rms} [mm]');
ylabel('y_{rms} [mm]');
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
print('thomx_errors_orbitrms.png','-dpng','-r300')

%%

dbetax_dipA=squeeze(DIP_align.E.dbeta(:,1,:));sdbetax_dipA=std(dbetax_dipA');maxdbetax_dipA=max(dbetax_dipA');dbetaxmax_dipA=max(abs(dbetax_dipA));dbetaxrms_dipA=std(dbetax_dipA);
dbetaz_dipA=squeeze(DIP_align.E.dbeta(:,2,:));sdbetaz_dipA=std(dbetaz_dipA');maxdbetaz_dipA=max(dbetaz_dipA');dbetazmax_dipA=max(abs(dbetaz_dipA));dbetazrms_dipA=std(dbetaz_dipA);
dbetax_dipF=squeeze(DIP_field.E.dbeta(:,1,:));sdbetax_dipF=std(dbetax_dipF');maxdbetax_dipF=max(dbetax_dipF');dbetaxmax_dipF=max(abs(dbetax_dipF));dbetaxrms_dipF=std(dbetax_dipF);
dbetaz_dipF=squeeze(DIP_field.E.dbeta(:,2,:));sdbetaz_dipF=std(dbetaz_dipF');maxdbetaz_dipF=max(dbetaz_dipF');dbetazmax_dipF=max(abs(dbetaz_dipF));dbetazrms_dipF=std(dbetaz_dipF);


dbetax_quadA=squeeze(QUAD_align.E.dbeta(:,1,:));sdbetax_quadA=std(dbetax_quadA');maxdbetax_quadA=max(dbetax_quadA');dbetaxmax_quadA=max(abs(dbetax_quadA));dbetaxrms_quadA=std(dbetax_quadA);
dbetaz_quadA=squeeze(QUAD_align.E.dbeta(:,2,:));sdbetaz_quadA=std(dbetaz_quadA');maxdbetaz_quadA=max(dbetaz_quadA');dbetazmax_quadA=max(abs(dbetaz_quadA));dbetazrms_quadA=std(dbetaz_quadA);
dbetax_quadF=squeeze(QUAD_field.E.dbeta(:,1,:));sdbetax_quadF=std(dbetax_quadF');maxdbetax_quadF=max(dbetax_quadF');dbetaxmax_quadF=max(abs(dbetax_quadF));dbetaxrms_quadF=std(dbetax_quadF);
dbetaz_quadF=squeeze(QUAD_field.E.dbeta(:,2,:));sdbetaz_quadF=std(dbetaz_quadF');maxdbetaz_quadF=max(dbetaz_quadF');dbetazmax_quadF=max(abs(dbetaz_quadF));dbetazrms_quadF=std(dbetaz_quadF);

dbetax_sextA=squeeze(SEXT_align.E.dbeta(:,1,:));sdbetax_sextA=std(dbetax_sextA');maxdbetax_sextA=max(dbetax_sextA');dbetaxmax_sextA=max(abs(dbetax_sextA));dbetaxrms_sextA=std(dbetax_sextA);
dbetaz_sextA=squeeze(SEXT_align.E.dbeta(:,2,:));sdbetaz_sextA=std(dbetaz_sextA');maxdbetaz_sextA=max(dbetaz_sextA');dbetazmax_sextA=max(abs(dbetaz_sextA));dbetazrms_sextA=std(dbetaz_sextA);
dbetax_sextF=squeeze(SEXT_field.E.dbeta(:,1,:));sdbetax_sextF=std(dbetax_sextF');maxdbetax_sextF=max(dbetax_sextF');dbetaxmax_sextF=max(abs(dbetax_sextF));dbetaxrms_sextF=std(dbetax_sextF);
dbetaz_sextF=squeeze(SEXT_field.E.dbeta(:,2,:));sdbetaz_sextF=std(dbetaz_sextF');maxdbetaz_sextF=max(dbetaz_sextF');dbetazmax_sextF=max(abs(dbetaz_sextF));dbetazrms_sextF=std(dbetaz_sextF);
    

figure(13);
plot(dbetaxmax_dipA,dbetazmax_dipA,'bd','MarkerSize',8,'DisplayName', 'Dipole alignment errors');
hold on
plot(dbetaxmax_dipF,dbetazmax_dipF,'rd','MarkerSize',8,'DisplayName','Dipole field errors');
plot(dbetaxmax_quadA,dbetazmax_quadA,'md','MarkerSize',8,'DisplayName', 'Quadrupole alignment errors');
plot(dbetaxmax_quadF,dbetazmax_quadF,'gd','MarkerSize',8,'DisplayName','Quadrupole field errors');
plot(dbetaxmax_sextA,dbetazmax_sextA,'cd','MarkerSize',8,'DisplayName', 'Sextupole alignment errors');
hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('[\Delta\beta_x/\beta_x]_{max}');
ylabel('[\Delta\beta_z/\beta_z]_{max}');
% xlim([-0.05 1.])
% ylim([-0.05 1])
xlim([-0.004 0.2])
ylim([-0.004 0.2])
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
print('thomx_errors_betabeatmax.png','-dpng','-r300')

figure(14);
plot(dbetaxrms_dipA,dbetazrms_dipA,'bo','MarkerSize',8,'DisplayName', 'Dipole alignment errors');
hold on
plot(dbetaxrms_dipF,dbetazrms_dipF,'ro','MarkerSize',8,'DisplayName','Dipole field errors');
plot(dbetaxrms_quadA,dbetazrms_quadA,'mo','MarkerSize',8,'DisplayName', 'Quadrupole alignment errors');
plot(dbetaxrms_quadF,dbetazrms_quadF,'go','MarkerSize',8,'DisplayName','Quadrupole field errors');
plot(dbetaxrms_sextA,dbetazrms_sextA,'co','MarkerSize',8,'DisplayName', 'Sextupole alignment errors');
hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('[\Delta\beta_x/\beta_x]_{rms}');
ylabel('[\Delta\beta_z/\beta_z]_{rms}');
% xlim([-0.01 0.35])
% ylim([-0.02 0.75])
xlim([-0.001 0.06])
ylim([-0.002 0.1])
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
print('thomx_errors_betabeatrms.png','-dpng','-r300')

%%

figure(15);
plot(DIP_align.E.SPos,sdbetax_dipA*100,'-r', 'Linewidth', 1.6);
hold on
plot(DIP_align.E.SPos,sdbetaz_dipA*100,'-b', 'Linewidth', 1.6);
plot(DIP_field.E.SPos,sdbetax_dipF*100,'--r', 'Linewidth', 1.6);
plot(DIP_field.E.SPos,sdbetaz_dipF*100,'--b', 'Linewidth', 1.6);
hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('s-position [m]');               
ylabel('[\Delta\beta/\beta]_{rms} [%]')
%legend('Horizontal','Vertical')
% u = legend('show','Location','NorthEast');
% set(u, 'FontSize',14)
u = legend([{'DIP alignment errors (Hor.)', 'DIP alignment errors (Vert.)', 'DIP field errors (Hor.)','DIP field errors (Vert.)'}]);
set(u,'Location','NorthEast','Orientation','vertical')
addlabel(1, 0, datestr(clock,0))
print('thomx_errors_betabeatrms_vs_s_DIP.png','-dpng','-r300')

figure(16);
plot(QUAD_align.E.SPos,sdbetax_quadA*100,'-r', 'Linewidth', 1.6);
hold on
plot(QUAD_align.E.SPos,sdbetaz_quadA*100,'-b', 'Linewidth', 1.6);
plot(QUAD_field.E.SPos,sdbetax_quadF*100,'--r', 'Linewidth', 1.6);
plot(QUAD_field.E.SPos,sdbetaz_quadF*100,'--b', 'Linewidth', 1.6);
hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('s-position [m]');               
ylabel('[\Delta\beta/\beta]_{rms} [%]')
% legend('Horizontal','Vertical')
% u = legend('show','Location','NorthEast');
% set(u, 'FontSize',14)
u = legend([{'QUAD alignment errors (Hor.)', 'QUAD alignment errors (Vert.)', 'QUAD field errors (Hor.)','QUAD field errors (Vert.)'}]);
set(u,'Location','NorthEast','Orientation','vertical')
addlabel(1, 0, datestr(clock,0))
print('thomx_errors_betabeatrms_vs_s_QUAD.png','-dpng','-r300')

figure(17);
plot(SEXT_align.E.SPos,sdbetax_sextA*100,'-r', 'Linewidth', 1.6);
hold on
plot(SEXT_align.E.SPos,sdbetaz_sextA*100,'-b', 'Linewidth', 1.6);
plot(SEXT_field.E.SPos,sdbetax_sextF*100,'--r', 'Linewidth', 1.6);
plot(SEXT_field.E.SPos,sdbetaz_sextF*100,'--b', 'Linewidth', 1.6);
hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('s-position [m]');               
ylabel('[\Delta\beta/\beta]_{rms} [%]')
% legend('Horizontal','Vertical')
% u = legend('show','Location','NorthEast');
% set(u, 'FontSize',14)
u = legend([{'SEXT alignment errors (Hor.)', 'SEXT alignment errors (Vert.)', 'SEXT field errors (Hor.)','SEXT field errors (Vert.)'}]);
set(u,'Location','NorthEast','Orientation','vertical')
addlabel(1, 0, datestr(clock,0))
print('thomx_errors_betabeatrms_vs_s_SEXT.png','-dpng','-r300')

%%


figure(18);
plot(DIP_align.E.SPos,maxdbetax_dipA*100,'-r', 'Linewidth', 1.6);
hold on
plot(DIP_align.E.SPos,maxdbetaz_dipA*100,'-b', 'Linewidth', 1.6);
plot(DIP_field.E.SPos,maxdbetax_dipF*100,'--r', 'Linewidth', 1.6);
plot(DIP_field.E.SPos,maxdbetaz_dipF*100,'--b', 'Linewidth', 1.6);
hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('s-position [m]');               
ylabel('[\Delta\beta/\beta]_{max} [%]')
u = legend([{'DIP alignment errors (Hor.)', 'DIP alignment errors (Vert.)', 'DIP field errors (Hor.)','DIP field errors (Vert.)'}]);
set(u,'Location','NorthEast','Orientation','vertical')
addlabel(1, 0, datestr(clock,0))
print('thomx_errors_betabeatmax_vs_s_DIP.png','-dpng','-r300')

figure(19);
plot(QUAD_align.E.SPos,maxdbetax_quadA*100,'-r', 'Linewidth', 1.6);
hold on
plot(QUAD_align.E.SPos,maxdbetaz_quadA*100,'-b', 'Linewidth', 1.6);
plot(QUAD_field.E.SPos,maxdbetax_quadF*100,'--r', 'Linewidth', 1.6);
plot(QUAD_field.E.SPos,maxdbetaz_quadF*100,'--b', 'Linewidth', 1.6);
hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('s-position [m]');               
ylabel('[\Delta\beta/\beta]_{max} [%]')
u = legend([{'QUAD alignment errors (Hor.)', 'QUAD alignment errors (Vert.)', 'QUAD field errors (Hor.)','QUAD field errors (Vert.)'}]);
set(u,'Location','NorthEast','Orientation','vertical')
addlabel(1, 0, datestr(clock,0))
print('thomx_errors_betabeatmax_vs_s_QUAD.png','-dpng','-r300')

figure(20);
plot(SEXT_align.E.SPos,maxdbetax_sextA*100,'-r', 'Linewidth', 1.6);
hold on
plot(SEXT_align.E.SPos,maxdbetaz_sextA*100,'-b', 'Linewidth', 1.6);
plot(SEXT_field.E.SPos,maxdbetax_sextF*100,'--r', 'Linewidth', 1.6);
plot(SEXT_field.E.SPos,maxdbetaz_sextF*100,'--b', 'Linewidth', 1.6);
hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('s-position [m]');               
ylabel('[\Delta\beta/\beta]_{max} [%]')
u = legend([{'SEXT alignment errors (Hor.)', 'SEXT alignment errors (Vert.)', 'SEXT field errors (Hor.)','SEXT field errors (Vert.)'}]);
set(u,'Location','NorthEast','Orientation','vertical')
addlabel(1, 0, datestr(clock,0))
print('thomx_errors_betabeatmax_vs_s_SEXT.png','-dpng','-r300')

%%

figure(21);
plot(DIP_align.E.tunes(:,1),DIP_align.E.tunes(:,2),'ob','MarkerSize',8,'DisplayName', 'Dipole alignment errors');
hold on
plot(DIP_field.E.tunes(:,1),DIP_field.E.tunes(:,2),'or','MarkerSize',8,'DisplayName', 'Dipole field errors');
plot(QUAD_align.E.tunes(:,1),QUAD_align.E.tunes(:,2),'om','MarkerSize',8,'DisplayName', 'Quadrupole alignment errors');
plot(QUAD_field.E.tunes(:,1),QUAD_field.E.tunes(:,2),'og','MarkerSize',8,'DisplayName', 'Quadrupole field errors');
plot(SEXT_align.E.tunes(:,1),SEXT_align.E.tunes(:,2),'oc','MarkerSize',8,'DisplayName', 'Sextupole alignment errors');
%plot(SEXT_field.E.tunes(:,1),SEXT_field.E.tunes(:,2),'xk','MarkerSize',8,'DisplayName', 'Sextupole field errors');
plot(DIP_align.E.tunes0(:,1),DIP_align.E.tunes0(:,2),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',9,'DisplayName', 'Initial tunes (3.17/1.64)')
hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
% u=legend('Initial','Trials')
% set(u, 'Location','SouthEast')
xlabel('\nu_x');
ylabel('\nu_z');
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
print('thomx_errors_tunes.png','-dpng','-r300')


%%

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(DIP_align.E.ringerrors(:,500),'comment',[],@plClosedOrbit)


figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(DIP_align.E.ringerrors(:,500),'comment',[],@pltmisalignments);


