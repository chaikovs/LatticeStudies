
clear all; close all; clc;

ring = ThomX_017_064_r56_02_chro00_AT2();
% ring_quadFF =atsetfieldvalues(ring,find(atgetcells(ring,'Class','Quadrupole')),...
%     'PassMethod','QuadMPoleFringePass' );
% ring = ring_quadFF;

%all_errors = load('data_ALLerrorsDA_align_err_100um_05mrad_field_err_0005');
%all_errors = load('data_ALLerrorsDA_align_err_30um_02mrad_field_err_0001');
% all_errors = load('data_ALLerrorsDA_align_err_30um_02mrad_field_err_0001_dpp_001');
% all_errors = load('data_ALLerrorsDA_align_err_30um_02mrad_field_err_0001_dpp_m001');
% all_errors = load('data_ALLerrorsDA_align_err_30um_02mrad_field_err_0001_dpp_0');
all_errors = load('data_ALLerrorsMultipBPM_align_err_30um_02mrad_field_err_0001_dpp_0');
% all_errors = load('data_ALLerrorsMultipBPM_align_err_30um_02mrad_field_err_0001_dpp_0');
%all_errors = load('data_ALLerrorsMultipBPM_align_err_100um_05mrad_field_err_0005_dpp_0_GOOD');

%all_errors = load('data_ALLerrorsMultipBPM_align_err_30um_02mrad_field_err_0001_dpp_0015');
%all_errors = load('data_ALLerrorsMultipBPM_align_err_30um_02mrad_field_err_0001_dpp_m0015');
%all_errors = load('data_ALLerrorsMultipBPM_align_err_30um_02mrad_field_err_0001_dpp_0_100seed');
rerr = all_errors.E.ringerrors(:,1);

%%
% get indices

%indBPM=find(atgetcells(ring,'Class','Monitor'));
indBPM=find(atgetcells(ring,'FamName','BPMx'));

indHCorT=find(atgetcells(ring,'Class','Corrector'));
indHCor = indHCorT(1:2:end);
indVCorT=find(atgetcells(ring,'Class','Corrector'));
indVCor = indVCorT(1:2:end);

sBPM=findspos(ring,indBPM);
sRING=findspos(ring,1:length(ring)+1);
%%

if ~exist('ModelRM')
ModelRM...
        =getresponsematrices_ira(...
        ring,...
        indBPM,...
        indHCor,...
        indVCor,...
        [],...
        [],...
        [],...
        0,...%[0 0 0 0 0 0]',...
        [1 2]); %      [1 2 3]);

end

%%
da_nturns = 100;
dpp = 0.0;

[lindata0, tunes0, chrom0] = twissring(ring, 0, 1:length(ring)+1,'chrom', 1e-8); % to get the tunes
beta0=cat(1,lindata0.beta);
SPos=cat(1,lindata0.SPos);

[xxda0,zzda0]=atdynap(ring,da_nturns,dpp,0.02);

%%
Nmachine = 500;
for kerr=1:Nmachine


% Get linear data with errors
    [lindataEr, nu, ch] = twissring(all_errors.E.ringerrors(:,kerr), 0, 1:length(all_errors.E.ringerrors(:,kerr))+1,'chrom', 1e-8); % to get the tunes
    betaEr  =cat(1,lindataEr.beta);
    tunesEr(kerr,:)=nu;
    chromEr(kerr,:)=ch;
    %
    dbetaEr(:,:,kerr)=(betaEr-beta0)./beta0;
    
    % DA
    try
        [xxerr,zzerr]=atdynap(all_errors.E.ringerrors(:,kerr),da_nturns,dpp,0.02);
    catch
        warning('Problem using function.  Assigning [xx,zz] value of 0.');
        xxerr = 0;
        zzerr = 0;
    end
    xxerrda(:,kerr)=xxerr;
    zzerrda(:,kerr)=zzerr;

% %DA
%     [xxerr,zzerr]=atdynap(all_errors.E.ringerrors(:,kerr),da_nturns,dpp,0.02);
%     xxerrda(:,kerr)=xxerr;
%     zzerrda(:,kerr)=zzerr;


[rcor,inCOD,hs,vs]=atcorrectorbit_ira(all_errors.E.ringerrors(:,kerr),...
    indBPM,...
    indHCor',...
    indVCor',...
    0,...%[0 0 0 0 0 0]',...
    [10 10],...
    [false true],...
    1.0,...
    ModelRM,...
    zeros(2,length(indBPM)),...
    [],...
    true);

[rcor,inCOD,hs,vs]=atcorrectorbit_ira(rcor,...
    indBPM,...
    indHCor',...
    indVCor',...
    0,...%[0 0 0 0 0 0]',...
    [10 10],...
    [false true],...
    1.0,...
    ModelRM,...
    zeros(2,length(indBPM)),...
    [],...
    true);

rcor_struct(:,kerr)  = rcor;

maxHS(kerr,:)= max(abs(hs));
maxVS(kerr,:)= max(abs(vs));

orbit6err = findorbit4Err(all_errors.E.ringerrors(:,kerr), inCOD, indBPM); % 1:length(rcor)+1
% xorbiterr{kerr}=orbit6err(1,:);
% yorbiterr{kerr}=orbit6err(3,:);
xorbiterr(:,kerr)=orbit6err(1,:);
yorbiterr(:,kerr)=orbit6err(3,:);

orbit6cor = findorbit4Err(rcor,inCOD,indBPM); %indBPM
xorbitcor(:,kerr)=orbit6cor(1,:);
yorbitcor(:,kerr)=orbit6cor(3,:);


% Get linear data
    [lindata, nu, ch] = twissring(rcor, 0, 1:length(rcor)+1,'chrom', 1e-8); % to get the tunes
    beta  =cat(1,lindata.beta);
    alpha = cat(1,lindata.alpha);
    disp  = cat(2,lindata.Dispersion);
    tunes(kerr,:)=nu;
    chrom(kerr,:)=ch;
    %
    dbeta(:,:,kerr)=(beta-beta0)./beta0;

 % DA
    [xx,zz]=atdynap(rcor,da_nturns,dpp,0.02);
    xxda(:,kerr)=xx;
    zzda(:,kerr)=zz;

end

%%

% nonzeroelem = sum(xxerrda~=0,1);
% problem_ind = find(nonzeroelem~=33);
% if ~isempty(problem_ind)
% % xxda=[xxda(:, 1:problem_ind - 1) xxda(:, problem_ind + 1:end)];
% % zzda=[zzda(:, 1:problem_ind - 1) zzda(:, problem_ind + 1:end)];
%   index = true(1, size(xxerrda, 2));
%   index([problem_ind]) = false;
%   xxerrda = xxerrda(:, index);
%   zzerrda = zzerrda(:, index);
% end


figure(111);
plot(xxda0*1e3,zzda0*1e3,'b-','LineWidth',2);hold on,
plot(xxerrda(:,1)*1e3,zzerrda(:,1)*1e3,'color',[0,0,0]+0.5,'LineWidth',1);
plot(xxda(:,1)*1e3,zzda(:,1)*1e3,'color',[0.1,0.5,0],'LineWidth',1);
plot(xxerrda*1e3,zzerrda*1e3,'color',[0,0,0]+0.5,'LineWidth',1);
plot(xxda*1e3,zzda*1e3,'color',[0.1,0.5,0],'LineWidth',1);
plot(xxda0*1e3,zzda0*1e3,'b-','LineWidth',2);hold off
legend('Initial: no errors','Trials with errors','Trials after orbit correction')
%hl = findobj(hobj,'type','line');
%set(hl,'LineWidth',3);
% ht = findobj(hobj,'type','text')
% set(ht,'FontSize',12);
% xlim([-50 50])
% ylim([0 30])
xlim([-30 30])
ylim([0 23])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x [mm]');                 % Add labels
ylabel('z [mm]');
title('DA (\deltap = 0%)')
%print('thomx_ALLminErrors_orbitCor_DA_init_seeds_dpp_0_multip100.png','-dpng','-r300')


%%

figure(112);
plot(tunes(:,1),tunes(:,2),'ob','MarkerSize',8);
hold on
plot(all_errors.E.tunes(:,1),all_errors.E.tunes(:,2),'go','MarkerSize',8);
plot(tunes0(:,1),tunes0(:,2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',9)
hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
u=legend('Tunes after orbit correction','Tunes with errors','Initial tunes (3.17/1.64)')
set(u, 'Location','NorthEast')
xlabel('\nu_x');
ylabel('\nu_z');
addlabel(1, 0, datestr(clock,0))

%%
% 
% % Stat over DA
% xxdai=[-60:0.2:60]*1e-3;
% zzdai0=interp1(xxda0,zzda0,xxdai','pchip',0);
% for kerr=1:nerr
%     zzdai(:,kerr)=interp1(xxda(:,kerr),zzda(:,kerr),xxdai','pchip',0);
% end
% 
% mzzdai=mean(zzdai,2);
% szzdai=std(zzdai')';
% DA_surf=sum(zzdai)*0.1e-3; % surface in m^2
% DA_surf0=sum(zzdai0)*0.1e-3; % surface in m^2
% 
% figure(112)
% plot(xxdai*1e3,zzdai0*1e3,'-b','LineWidth',2);hold on
% plot(xxdai*1e3,mzzdai*1e3,'-r','LineWidth',2);
% plot(xxdai*1e3,(mzzdai-szzdai)*1e3,'--m','LineWidth',2);
% plot(xxdai*1e3,(mzzdai+szzdai)*1e3,'--m','LineWidth',2);hold off
% legend('Initial: no errors','Mean of 100 seeds','Mean - \sigma','Mean + \sigma')
% % xlim([-60 60])
% % ylim([0 25])
% xlim([-30 30])
% ylim([0 20])
% grid on
% set(gcf,'color','w')
% set(gca,'fontsize',16');
% xlabel('x [mm]');                 % Add labels
% ylabel('z [mm]');
% title('DA (\deltap = 1%)')
% %print('thomx_ALLminErrors_DA_init_,meanseeds1_dpp_001.png','-dpng','-r300')
% 
% figure(113)
% plot(xxda0*1e3,zzda0*1e3,'b-','LineWidth',2);hold on
% plot(xxdai*1e3,mzzdai*1e3,'-r','LineWidth',2);
% plot(xxda*1e3,zzda*1e3,'color',[0,0,0]+0.5,'LineWidth',1);
% plot(xxdai*1e3,mzzdai*1e3,'-r','LineWidth',2);
% plot(xxda0*1e3,zzda0*1e3,'b-','LineWidth',2)
% % plot(xxda0*1e3,zzda0*1e3,'m-','LineWidth',2)
% %plot(xxdai*1e3,(mzzdai-szzdai)*1e3,'--m','LineWidth',2);
% %plot(xxdai*1e3,(mzzdai+szzdai)*1e3,'--m','LineWidth',2);hold off
% legend('Initial: no errors','Mean of 100 seeds')
% xlim([-60 60])
% ylim([0 30])
% grid on
% set(gcf,'color','w')
% set(gca,'fontsize',16');
% xlabel('x [mm]');                 % Add labels
% ylabel('z [mm]');
% title('DA (\deltap = 0)')
% %print('thomx_ALLminErrors_DA_init_meanseeds2.png','-dpng','-r300')
% 
% figure(114);
% plot(DA_surf0*1e6, 0 ,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10);hold on
% h=histogram(DA_surf*1e6); 
% h.FaceColor = [0 0.5 0.5];
% hold off    
% set(gcf,'color','w')
% set(gca,'fontsize',16');
% xlabel('DA surface [mm^2]');
% ylabel('Entries');
% title('DA (\deltap = 1%)')
% u=legend('Initial: no errors','Trials with errors')
% set(u,'Location','NorthWest','Orientation','vertical')
% %print('thomx_ALLminErrors_DA_surf_dpp_001.png','-dpng','-r300')

%%

%   
% [rcor,inCOD,hs,vs]=atcorrectorbit(rerr,...
%     indBPM,...
%     indHCor',...
%     indVCor',...
%     [0 0 0 0 0 0]',...
%     [12 12],...
%     [false true],...
%     1.0,...
%     ModelRM,...
%     zeros(2,length(indBPM)),...
%     [],...
%     true);

%%

% %find new closed orbit
% o=findorbit6Err(rerr,indBPM,inCOD);
% oxe=o(1,:);
% oye=o(3,:);
% 
% o=findorbit6Err(rcor,indBPM,inCOD);
% oxc=o(1,:);
% oyc=o(3,:);
% 


% figure(11);
% subplot(2,1,1);
% plot(sBPM,oxe,'.-');hold on; plot(sBPM,oxc,'.-');
% legend('before','after');
% xlabel('s [m]');
% ylabel('hor. COD');
% subplot(2,1,2);
% plot(sBPM,oye,'.-');hold on; plot(sBPM,oyc,'.-');
% legend('before','after');
% xlabel('s [m]');
% ylabel('ver. COD');
% %saveas(gca,'OrbitCor.fig');
% %export_fig('OrbitCor.jpg','-r300');

% xorbiterr_plotT=[xorbiterr{:}];
% xorbiterr_plot = reshape(xorbiterr_plotT, [12, Nmachine]);
% yorbiterr_plotT=[yorbiterr{:}];
% yorbiterr_plot = reshape(yorbiterr_plotT, [12, Nmachine]);
% xorbitcor_plotT = [xorbitcor{:}];
% xorbitcor_plot = reshape(xorbitcor_plotT, [12, Nmachine]);
% yorbitcor_plotT = [yorbitcor{:}];
% yorbitcor_plot = reshape(yorbitcor_plotT, [12, Nmachine]);
% 

figure(12);
plot(sBPM,1e3.*xorbiterr(:,1:50),'b.-','MarkerSize',10)
hold on; 
%plot(sBPM,1e3.*xorbiterr(indBPM,:),'kx','MarkerSize',15)
plot(sBPM,1e3.*xorbitcor(:,1:50),'r.-');
%plot(sBPM,1e3.*xorbitcor(indBPM,:),'kx','MarkerSize',15);
set(gcf,'color','w')
set(gca,'fontsize',16');
 xlim([0 18])
% ylim([0 6])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('s-position [m]');
ylabel('x [mm]');
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLminErrors_orbitCor_x_RING_multip.png','-dpng','-r300')

% figure(121);
% plot(sBPM,1e3.*xorbiterr(indBPM,:),'bd-','MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','k')
% hold on; 
% plot(sBPM,1e3.*xorbitcor(indBPM,:),'rd-','MarkerSize',6);
% set(gcf,'color','w')
% set(gca,'fontsize',16');
%  xlim([0 18])
% % ylim([0 6])
% grid on
% set(gcf,'color','w')
% set(gca,'fontsize',16');
% xlabel('s-position [m]');
% ylabel('x [mm]');
% addlabel(1, 0, datestr(clock,0))
% %print('thomx_ALLmaxErrors_orbitCor_x_BPM.png','-dpng','-r300')

figure(13);
plot(sBPM,1e3*yorbiterr(:,1:50),'b.-','MarkerSize',8)
hold on; 
%plot(sBPM,1e3*yorbiterr(indBPM,:),'b.-','MarkerSize',8)
plot(sBPM,1e3*yorbitcor(:,1:50),'r.-');
%plot(sBPM,1e3*yorbitcor(indBPM,:),'r.-');
set(gcf,'color','w')
set(gca,'fontsize',16');
 xlim([0 18])
% ylim([0 6])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('s-position [m]');
ylabel('y [mm]');
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLminErrors_orbitCor_y_RING_multip.png','-dpng','-r300')

% figure(131);
% plot(sBPM,1e3*yorbiterr(indBPM,:),'bd-','MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','k')
% hold on; 
% plot(sBPM,1e3*yorbitcor(indBPM,:),'rd-','MarkerSize',6);
% set(gcf,'color','w')
% set(gca,'fontsize',16');
%  xlim([0 18])
% % ylim([0 6])
% grid on
% set(gcf,'color','w')
% set(gca,'fontsize',16');
% xlabel('s-position [m]');
% ylabel('y [mm]');
% addlabel(1, 0, datestr(clock,0))
% %print('thomx_ALLmaxErrors_orbitCor_y_BPM.png','-dpng','-r300')

%%

figure(132);
plot(1e3*max(abs(xorbiterr)),1e3*max(abs(yorbiterr)),'+b','MarkerSize',8,'DisplayName', 'Before orbit correction')
hold on
plot(1e3*max(abs(xorbitcor)),1e3*max(abs(yorbitcor)),'+r','MarkerSize',8,'DisplayName', 'After orbit correction')
hold off
set(gcf,'color','w')
set(gca,'fontsize',16');
% xlim([-0.5 25])
% ylim([-0.5 15])
xlim([-0.1 8])
ylim([-0.1 5])
grid on
xlabel('x_{max} [mm]');
ylabel('y_{max} [mm]');
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLminErrors_orbitCor_max_multip.png','-dpng','-r300')

figure(133);
plot(1e3*std(xorbiterr),1e3*std(yorbiterr),'+b','MarkerSize',8,'DisplayName', 'Before orbit correction')
hold on
plot(1e3*std(xorbitcor),1e3*std(yorbitcor),'+r','MarkerSize',8,'DisplayName', 'After orbit correction')
hold off
set(gcf,'color','w')
set(gca,'fontsize',16');
% xlim([-0.5 15])
% ylim([-0.5 10])
xlim([-0.05 4])
ylim([-0.05 3])
grid on
xlabel('x_{rms} [mm]');
ylabel('y_{rms} [mm]');
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLminErrors_orbitCor_rms_multip.png','-dpng','-r300')

%%

Lh=atgetfieldvalues(rcor,indHCor,'Length');
Lv=atgetfieldvalues(rcor,indVCor,'Length');
hsL=maxHS.*Lh(1);
vsL=maxVS.*Lv(1);

figure(14);
histogram(hsL*1e3,10)
set(gcf,'color','w')
set(gca,'fontsize',16');
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('Maximal horizontal kick [mrad]');
ylabel('Number of machines');
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLminErrors_orbitCor_HCorStrength_multip.png','-dpng','-r300')


figure(15);
histogram(vsL*1e3,10)
set(gcf,'color','w')
set(gca,'fontsize',16');
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('Maximal vertical kick [mrad]');
ylabel('Number of machines');
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLminErrors_orbitCor_VCorStrength_multip.png','-dpng','-r300')


% figure(141);
% histogram(hsL*1e3)
% set(gcf,'color','w')
% set(gca,'fontsize',16');
% grid on
% set(gcf,'color','w')
% set(gca,'fontsize',16');
% xlabel('Horizontal kick [mrad]');
% ylabel('Entries');
% addlabel(1, 0, datestr(clock,0))
% 
% figure(142);
% histogram(vsL*1e3)
% set(gcf,'color','w')
% set(gca,'fontsize',16');
% grid on
% set(gcf,'color','w')
% set(gca,'fontsize',16');
% xlabel('Vertical kick [mrad]');
% ylabel('Entries');
% addlabel(1, 0, datestr(clock,0))
% 

% plot output
% 
% figure(15);
% subplot(2,1,1);bar(hs);ylabel('hor.')
% subplot(2,1,2);bar(vs);ylabel('ver.')

%%

dbetax=squeeze(dbeta(:,1,:));sdbetax=std(dbetax');maxdbetax=max(dbetax');dbetaxmax=max(abs(dbetax));dbetaxrms=std(dbetax);
dbetaz=squeeze(dbeta(:,2,:));sdbetaz=std(dbetaz');maxdbetaz=max(dbetaz');dbetazmax=max(abs(dbetaz));dbetazrms=std(dbetaz);

dbetaxEr=squeeze(dbetaEr(:,1,:));sdbetaxEr=std(dbetaxEr');maxdbetaxEr=max(dbetaxEr');dbetaxmaxEr=max(abs(dbetaxEr));dbetaxrmsEr=std(dbetaxEr);
dbetazEr=squeeze(dbetaEr(:,2,:));sdbetazEr=std(dbetazEr');maxdbetazEr=max(dbetazEr');dbetazmaxEr=max(abs(dbetazEr));dbetazrmsEr=std(dbetazEr);

figure(16);
plot(SPos,sdbetax*100,'-r', 'Linewidth', 1.6);hold on
plot(SPos,sdbetaz*100,'-b', 'Linewidth', 1.6);hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
legend('Horizontal','Vertical')
xlabel('s-position [m]');               
ylabel('[\Delta\beta/\beta]_{rms} [%]')
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLminErrors_orbitCor_betabeatrms_vs_s_multip.png','-dpng','-r300')

figure(17);
plot(SPos,maxdbetax*100,'-r', 'Linewidth', 1.6);hold on
plot(SPos,maxdbetaz*100,'-b', 'Linewidth', 1.6);hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
legend('Horizontal','Vertical')
xlabel('s-position [m]');               
ylabel('[\Delta\beta/\beta]_{max} [%]')
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLminErrors_orbitCor_betabeatmax_vs_s_multip.png','-dpng','-r300')

figure(18);
plot(dbetaxmax,dbetazmax,'rd','DisplayName', 'After orbit correction');
hold on
plot(dbetaxmaxEr,dbetazmaxEr,'bd','DisplayName', 'Before orbit correction');
hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('[\Delta\beta_x/\beta_x]_{max}');
ylabel('[\Delta\beta_z/\beta_z]_{max}');
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLminErrors_orbitCor_betabeatmax_multip.png','-dpng','-r300')

figure(19);
plot(dbetaxrms,dbetazrms,'ro','DisplayName', 'After orbit correction');
hold on
plot(dbetaxrmsEr,dbetazrmsEr,'bo','DisplayName', 'Before orbit correction');
hold off
axis tight
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('[\Delta\beta_x/\beta_x]_{rms}');
ylabel('[\Delta\beta_z/\beta_z]_{rms}');
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLminErrors_orbitCor_betabeatrms_multip.png','-dpng','-r300')

