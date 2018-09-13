
clear all; close all; clc;

%ring = ThomX_017_064_r56_02_chro00_AT2();

all_errors = load('data_ALLerrorsMultipBPM_align_err_100um_05mrad_field_err_0005_dpp_0_GOOD');

ring = ThomX_017_064_r56_02_chro00_multip_AT2();

%%


indm=find(atgetcells(ring,'FamName','BPMx'));
sBPM=findspos(ring,indm);

[lindata0, tunes0, chrom0] = twissring(ring, 0, 1:length(ring)+1,'chrom', 1e-8); % to get the tunes
beta0=cat(1,lindata0.beta);
SPos=cat(1,lindata0.SPos);
[xxda0,zzda0]=atdynap(ring,100,0.0,0.02);
sizebeta=size(beta0);
%%

figure(11);
plot(1e3*all_errors.E.xorbitmax,1e3*all_errors.E.yorbitmax,'+b','MarkerSize',9,'DisplayName','Without BPM errors')
hold on
plot(1e3*all_errors.E.xorbitmaxErr,1e3*all_errors.E.yorbitmaxErr,'+r','MarkerSize',9,'DisplayName','With BPM errors')
hold off
set(gcf,'color','w')
set(gca,'fontsize',16');
% xlim([0 5])
% ylim([0 5])
xlim([0 17])
ylim([0 13])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x_{max} [mm]');
ylabel('y_{max} [mm]');
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
print('thomx_ALLmaxErrors_orbitmax_dpp_0_multip.png','-dpng','-r300')


figure(12);
plot(1e3*all_errors.E.xorbitrms,1e3*all_errors.E.yorbitrms,'xb','MarkerSize',8,'DisplayName','Without BPM errors')
hold on
plot(1e3*all_errors.E.xorbitrmsErr,1e3*all_errors.E.yorbitrmsErr,'+r','MarkerSize',8,'DisplayName','With BPM errors')
hold off
set(gcf,'color','w')
set(gca,'fontsize',16');
% xlim([0 4])
% ylim([0 3])
xlim([0 12])
ylim([0 9])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x_{rms} [mm]');
ylabel('y_{rms} [mm]');
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
print('thomx_ALLmaxErrors_orbitrms_dpp_0_multip.png','-dpng','-r300')


figure(121);
plot(sBPM,1e3.*all_errors.E.xorbit(:, 1),'b.-','MarkerSize',10,'DisplayName','Orbit')
hold on; 
plot(sBPM,1e3.*all_errors.E.xorbitErr(:, 1),'rx-','MarkerSize',10,'DisplayName','BPM readings');
plot(sBPM,1e3.*all_errors.E.xorbit(:, 1:3),'b.-','MarkerSize',10,'HandleVisibility','off')
plot(sBPM,1e3.*all_errors.E.xorbitErr(:, 1:3),'rx-','MarkerSize',10,'HandleVisibility','off');
hold off
set(gcf,'color','w')
set(gca,'fontsize',16');
 xlim([0 18])
% ylim([0 6])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('s-position [m]');
ylabel('x [mm]');
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
print('thomx_ALLmaxErrors_xorbit_dpp_0_multip.png','-dpng','-r300')

figure(122);
plot(sBPM,1e3.*all_errors.E.yorbit(:, 1),'b.-','MarkerSize',10,'DisplayName','Orbit')
hold on; 
plot(sBPM,1e3.*all_errors.E.yorbitErr(:, 1),'rx-','MarkerSize',10,'DisplayName','BPM readings');
plot(sBPM,1e3.*all_errors.E.yorbit(:, 1:3),'b.-','MarkerSize',10,'HandleVisibility','off')
plot(sBPM,1e3.*all_errors.E.yorbitErr(:, 1:3),'rx-','MarkerSize',10,'HandleVisibility','off');
hold off
set(gcf,'color','w')
set(gca,'fontsize',16');
 xlim([0 18])
% ylim([0 6])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('s-position [m]');
ylabel('y [mm]');
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
print('thomx_ALLmaxErrors_yorbit_dpp_0_multip.png','-dpng','-r300')



dbetax=squeeze(all_errors.E.dbeta(:,1,:));sdbetax=std(dbetax');maxdbetax=max(dbetax');dbetaxmax=max(abs(dbetax));
dbetaz=squeeze(all_errors.E.dbeta(:,2,:));sdbetaz=std(dbetaz');maxdbetaz=max(dbetaz');dbetazmax=max(abs(dbetaz));

figure(13);
plot(SPos,sdbetax*100,'-r', 'Linewidth', 1.6);hold on
plot(SPos,sdbetaz*100,'-b', 'Linewidth', 1.6);hold off
axis tight
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
legend('Horizontal','Vertical')
xlabel('s-position [m]');               
ylabel('[\Delta\beta/\beta]_{rms} [%]')
addlabel(1, 0, datestr(clock,0))
print('thomx_ALLmaxErrors_betabeatrms_vs_s_multip.png','-dpng','-r300')

figure(14);
plot(SPos,maxdbetax*100,'-r', 'Linewidth', 1.6);hold on
plot(SPos,maxdbetaz*100,'-b', 'Linewidth', 1.6);hold off
axis tight
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
legend('Horizontal','Vertical')
xlabel('s-position [m]');               
ylabel('[\Delta\beta/\beta]_{max} [%]')
addlabel(1, 0, datestr(clock,0))
print('thomx_ALLmaxErrors_betabeatmax_vs_s_multip.png','-dpng','-r300')

figure(15);
plot(dbetaxmax,dbetazmax,'ro');
axis tight
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('[\Delta\beta_x/\beta_x]_{max}');
ylabel('[\Delta\beta_z/\beta_z]_{max}');
addlabel(1, 0, datestr(clock,0))
print('thomx_ALLmaxErrors_betabeatmax_dpp_0_multip.png','-dpng','-r300')

figure(16);
plot(all_errors.E.tunes(:,1),all_errors.E.tunes(:,2),'ob','MarkerSize',8);
hold on
plot(all_errors.E.tunes0(:,1),all_errors.E.tunes0(:,2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',9)
hold off
axis tight
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
u=legend('Random seeds','Initial tunes (3.17/1.64)')
set(u, 'Location','NorthEast')
xlabel('\nu_x');
ylabel('\nu_z');
addlabel(1, 0, datestr(clock,0))
print('thomx_ALLmaxErrors_tunes_dpp_0_multip.png','-dpng','-r300')


%% DA


figure(17);
plot(all_errors.E.xda0*1e3,all_errors.E.zda0*1e3,'b-','LineWidth',2);hold on,
plot(all_errors.E.xda*1e3,all_errors.E.zda*1e3,'color',[0,0,0]+0.5,'LineWidth',1);
plot(all_errors.E.xda0*1e3,all_errors.E.zda0*1e3,'b-','LineWidth',2);hold off
legend('Initial: no errors','Trials with errors')
% xlim([-50 50])
% ylim([0 30])
xlim([-30 30])
ylim([0 15])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x [mm]');                 % Add labels
ylabel('z [mm]');
title('DA (\deltap = 0%)')
print('thomx_ALLmaxErrors_DA_init_seeds_dpp_0_multip.png','-dpng','-r300')

% Stat over DA
xdai=[-40:0.2:40]*1e-3;
zdai0=interp1(all_errors.E.xda0,all_errors.E.zda0,xdai','pchip',0);
for kerr=1:length(all_errors.E.xda)%nerr
    zdai(:,kerr)=interp1(all_errors.E.xda(:,kerr),all_errors.E.zda(:,kerr),xdai','pchip',0);
end

mzdai=mean(zdai,2);
szdai=std(zdai')';
DA_surf=sum(zdai)*0.1e-3; % surface in m^2
DA_surf0=sum(zdai0)*0.1e-3; % surface in m^2

figure(18)
plot(xdai*1e3,zdai0*1e3,'-b','LineWidth',2);hold on
plot(xdai*1e3,mzdai*1e3,'-r','LineWidth',2);
plot(xdai*1e3,(mzdai-szdai)*1e3,'--m','LineWidth',2);
plot(xdai*1e3,(mzdai+szdai)*1e3,'--m','LineWidth',2);hold off
legend('Initial: no errors','Mean of 500 seeds','Mean - \sigma','Mean + \sigma')
% xlim([-60 60])
% ylim([0 25])
xlim([-30 30])
ylim([0 15])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x [mm]');                 % Add labels
ylabel('z [mm]');
title('DA (\deltap = 0%)')
print('thomx_ALLmaxErrors_DA_init_,meanseeds1_dpp_0_multip.png','-dpng','-r300')

figure(181)
plot(all_errors.E.xda0*1e3,all_errors.E.zda0*1e3,'b-','LineWidth',2);hold on
plot(xdai*1e3,mzdai*1e3,'-r','LineWidth',2);
plot(all_errors.E.xda*1e3,all_errors.E.zda*1e3,'color',[0,0,0]+0.5,'LineWidth',1);
plot(xdai*1e3,mzdai*1e3,'-r','LineWidth',2);
plot(all_errors.E.xda0*1e3,all_errors.E.zda0*1e3,'b-','LineWidth',2)
% plot(xxda0*1e3,zzda0*1e3,'m-','LineWidth',2)
%plot(xxdai*1e3,(mzzdai-szzdai)*1e3,'--m','LineWidth',2);
%plot(xxdai*1e3,(mzzdai+szzdai)*1e3,'--m','LineWidth',2);hold off
legend('Initial: no errors','Mean of 500 seeds')
xlim([-30 30])
ylim([0 15])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x [mm]');                 % Add labels
ylabel('z [mm]');
title('DA (\deltap = 0)')
print('thomx_ALLmaxErrors_DA_init_meanseeds2_multip.png','-dpng','-r300')

figure(19);
plot(DA_surf0*1e6, 0 ,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10);hold on
h=histogram(DA_surf*1e6); 
h.FaceColor = [0 0.5 0.5];
hold off    
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('DA surface [mm^2]');
ylabel('Entries');
title('DA (\deltap = 0%)')
u=legend('Initial: no errors','Trials with errors')
set(u,'Location','NorthEast','Orientation','vertical')
print('thomx_ALLmaxErrors_DA_surf_dpp_0_multip.png','-dpng','-r300')


%%

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(all_errors.E.ringerrors(:,500),'comment',[],@plClosedOrbit)


figure('units','normalized','position',[0.1 0.4 0.45 0.35])
atplot(all_errors.E.ringerrors(:,500),'comment',[],@pltmisalignments);
print('thomx_ALLmaxErrors_misalignments_dpp_0_multip.png','-dpng','-r300')
%export_fig('thomx_ALLminErrors_misalignments_dpp_0_multip.pdf','-r300');

return
