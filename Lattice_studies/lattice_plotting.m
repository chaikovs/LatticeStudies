global GLOBVAL 
GLOBVAL.E0=50E6;
GLOBVAL.LatticeFile='test';

dir ='/Users/ichaikov/Documents/MATLAB/thomx-mml/machine/THOMX/StorageRing/Lattices/';
RING = ThomX_017_064_r56_02_chro00;
%RING = ThomX_016_058_r56_02_chro22(); 
 
 %% Betatron phase

dpp = 0;
[lindata,tune,chrom]=atlinopt(RING,dpp,1:length(RING)+1); 
dispersion=cat(2,lindata.Dispersion)';
beta=cat(1,lindata.beta);
mu=cat(1,lindata.mu);
spos=cat(1,lindata.SPos);
%%
indm=find(atgetcells(RING,'FamName','BPMx'))
indqp1=find(atgetcells(RING,'FamName','QP1'))

 figure
 set(gcf,'color','w')
 plot(spos, mu(:,1)/(2*pi),'b-','LineWidth',2,'DisplayName', 'Horizontal betatron phase')
 hold on
 plot(spos(indm), mu(indm,1)/(2*pi),'ko','LineWidth',2,'DisplayName', 'Horizontal betatron phase')
 plot(spos(indqp1), mu(indqp1,1)/(2*pi),'r*','LineWidth',2,'DisplayName', 'Horizontal betatron phase')
  plot(spos, mu(:,2)/(2*pi),'r-','LineWidth',2,'DisplayName', 'Vertical betatron phase')
 hold off
 grid on
set(gca,'fontsize',20);
xlim([0 max(spos)])
xlabel('s-position [m]');                 % Add labels
ylabel('Betatron phase \mu');
u = legend('show','Location','NorthWest');
set(u,'FontSize',14)
%print('thomx_betatron_phase.png','-dpng','-r300')

%%

figure
set(gcf,'color','w')
h1 = subplot(5,1,[1 4]);
plot(spos, mu(:,1)/(2*pi),'b-','LineWidth',2,'DisplayName', 'Horizontal betatron phase')
hold on
plot(spos, mu(:,2)/(2*pi),'r-','LineWidth',2,'DisplayName', 'Vertical betatron phase')
hold off
grid on
set(gca,'fontsize',18);
xlim([0 max(spos)])
%set(gca,'xtick',[])
set(gca,'xticklabel',[])
%xlabel('s-position [m]');                 % Add labels
ylabel('Betatron phase \mu');
u = legend('show','Location','NorthWest');
set(u,'FontSize',14)
h2 = subplot(5,1,5);
drawlattice
set(h2,'YTick',[])
set(gca,'FontSize',18)
xlabel('s - position [m]');
linkaxes([h1 h2],'x')
set([h1 h2],'XGrid','On','YGrid','On');
print('thomx_betatron_phase_lat.png','-dpng','-r300')

%% Optics 1/4

figure(1)
h1 = subplot(5,1,[1 4]);
set(gca,'FontSize',18)
plot(spos,beta(:,1),'.-b', 'Markersize',10, 'Linewidth', 1.6);
hold on
plot(spos,beta(:,2),'.-r', 'Markersize',10, 'Linewidth', 1.6);
plot(spos, 10*dispersion(:,1),'.-g', 'Markersize',10, 'Linewidth', 1.6)
hold off
set(gca,'fontsize',18);
% xlim([0 S(end)]);
xlim([0 4.5]);
set(gca,'xticklabel',[])
ylabel('\beta_x,\beta_z,\eta_x [m]');
%title('Optical-functions');
u = legend({'\beta_x','\beta_z','10*\eta_x'});
set(u,'Location','NorthWest')
h2 = subplot(5,1,5);
set(gca,'FontSize',16)
drawlattice 
set(h2,'YTick',[])
set(gca,'FontSize',18)
xlabel('s - position [m]');
linkaxes([h1 h2],'x')
xlim([0 4.5]);
set([h1 h2],'XGrid','On','YGrid','On');
print('thomx_optics_period_WP2.png','-dpng','-r300')

%%

atplot(RING,[0 4.5])
print('thomx_optics_period_atplot_WP2.png','-dpng','-r300')

%%

atplot(RING)
print('thomx_optics_atplot_WP0.png','-dpng','-r300')
%%

dpp = 0.01;
[lindata_OffM,tune_OffM,chrom_OffM]=atlinopt(RING,dpp,1:length(RING)+1); 
dispersion_OffM=cat(2,lindata_OffM.Dispersion)';
beta_OffM=cat(1,lindata_OffM.beta);
spos_OffM=cat(1,lindata_OffM.SPos);

dpp = -0.01;
[lindata_OffMneg,tune_OffMneg,chrom_OffMneg]=atlinopt(RING,dpp,1:length(RING)+1); 
dispersion_OffMneg=cat(2,lindata_OffMneg.Dispersion)';
beta_OffMneg=cat(1,lindata_OffMneg.beta);
spos_OffMneg=cat(1,lindata_OffMneg.SPos);


%%

figure('units','normalized','position',[0.3 0.3 0.55 0.45])
h1 = subplot(5,1,[1 4]);
set(gca,'FontSize',18)
plot(spos,beta(:,1),'-k', 'Markersize',10, 'Linewidth', 1.6);
hold on
plot(spos,beta(:,2),'-k', 'Markersize',10, 'Linewidth', 1.6);
%plot(spos, 10*dispersion(:,1),'.-g', 'Markersize',10, 'Linewidth', 1.6)
plot(spos_OffMneg,beta_OffMneg(:,1),'--b', 'Markersize',10, 'Linewidth', 1.6);
plot(spos_OffMneg,beta_OffMneg(:,2),'--r', 'Markersize',10, 'Linewidth', 1.6);
plot(spos_OffM,beta_OffM(:,1),'-b', 'Markersize',10, 'Linewidth', 1.6);
plot(spos_OffM,beta_OffM(:,2),'-r', 'Markersize',10, 'Linewidth', 1.6);
%plot(spos_OffM, 10*dispersion_OffM(:,1),'.--g', 'Markersize',10, 'Linewidth', 1.6)
hold off
set(gca,'fontsize',18);
ylabel('\beta_x,\beta_z [m]');
u = legend({'\beta_x', '\beta_z', '\beta_x (dp/p = -1%)','\beta_z (dp/p = -1%)','\beta_x (dp/p = 1%)','\beta_z (dp/p = 1%)'});
set(u,'Location','NorthWest','Orientation','horizontal')
% xlim([0 S(end)]);
xlim([0 4.5]);
ylim([0 10]);
set(gca,'xticklabel',[])
h2 = subplot(5,1,5);
set(gca,'FontSize',16)
drawlattice 
set(h2,'YTick',[])
set(gca,'FontSize',18)
xlabel('s - position [m]');
linkaxes([h1 h2],'x')
xlim([0 4.5]);
set([h1 h2],'XGrid','On','YGrid','On');

print('thomx_optics_offmomentum001_WP0.png','-dpng','-r300')

%%
%  global THERING
%  [Dx, Dy, Sx, Sy] = modeldisp
% 
% 
% 
% [TD, tune] = twissring(THERING,0,1:(length(THERING)+1));
% BETA = cat(1,TD.beta);
% MU   = cat(1,TD.mu); % not normalized to 2pi
% S  = cat(1,TD.SPos);
% %disp(tune)
% 
% % Figure to check
% % plot betax and betay in two subplots
% figure(1)
% h1 = subplot(5,1,[1 4]);
% set(gca,'FontSize',18)
% plot(S,BETA(:,1),'.-b', 'Markersize',10, 'Linewidth', 1.6);
% hold on
% plot(S,BETA(:,2),'.-r', 'Markersize',10, 'Linewidth', 1.6);
% plot(Sx, 10*Dx,'.-g', 'Markersize',10, 'Linewidth', 1.6)
% hold off
% % xlim([0 S(end)]);
% xlim([0 4.5]);
% ylabel('\beta_x,\beta_z,\eta_x [m]');
% %title('Optical-functions');
% u = legend({'\beta_x','\beta_z','10*\eta_x'});
% set(u,'Location','NorthEast')
% h2 = subplot(5,1,5);
% set(gca,'FontSize',16)
% drawlattice 
% set(h2,'YTick',[])
% xlabel('s - position [m]');
% linkaxes([h1 h2],'x')
% xlim([0 4.5]);
%set([h1 h2],'XGrid','On','YGrid','On');

