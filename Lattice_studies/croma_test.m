% 
    global GLOBVAL
    GLOBVAL.E0=50E6;
    GLOBVAL.LatticeFile='test';

dir ='/Users/ichaikov/Documents/MATLAB/thomx-mml/machine/THOMX/StorageRing/Lattices/';
RING = ThomX_017_064_r56_02_chro00; 
%RING = ThomX_017_064_r56_02_chro00_multip_AT2();
%RING = ThomX_016_058_r56_02_chro22;
%RING = ThomX_017_064_r56_02_chro00_AT2;


% %QUAD FF
% ring_quadFF =atsetfieldvalues(RING,find(atgetcells(RING,'Class','Quadrupole')),...
%     'PassMethod','QuadMPoleFringePass' );
% RING = ring_quadFF;
% 
% RING=scalesext(RING,'SX1',0.);
% RING=fitchrom_alex(RING,[0.0 -0.0],'SX2' ,'SX3' );

atplot(RING)
%%
%L =   findspos(RING,length(RING)+1)
% atplot(RING)
% print('thomx_lattice.png','-dpng','-r300')

% dnu=[0.0 0.2];
% for k=1:5
% RING=fittune2_alex(RING, [(3.17+k*dnu(1)/5) (1.64+k*dnu(2)/5)], 'QP3', 'QP41'); 
% end

%RING=fittune2_alex(RING, [3.16 1.7], 'QP31', 'QP4'); % was 1.58
%RING=scalesext(RING,'SX1',1);
%RING=fitchrom_alex(RING,[2 -2],'SX2' ,'SX3' );

%RING=scalesext(RING,'SX1',0.);
%RING=fitchrom_alex(RING,[2.0 -2.0],'SX2' ,'SX3' );
%ind=find(atgetcells(RING,'Class','Quadrupole'));RING=atsetfieldvalues(RING,ind,'PassMethod','QuadMPoleFringePass' );
% RING=fittune2_alex(RING, [3.12 1.58], 'QP31', 'QP4');
%RING=fittune2_alex(RING, [3.16 1.58], 'QP31', 'QP4');


%atplot(RING)
%atplot_alex(RING)
%figure(101);atplot(RING,@plotRDT,'geometric1');

% rin=[1e-4; 0.00; 3e-3; 0; 0; 0];
% [X]=ringpass(RING,rin,100); 
% figure(50);plot(X(1,:),X(2,:),'.b')
% figure(51);plot(X(3,:),X(4,:),'.b')
% findtune_multi(reshape(X(1,:),1,[])',3)
% findtune_multi(reshape(X(3,:),1,[])',3)
%return

%%
% [xx,zz]=atdynap(RING,500,-0.0,0.02);
% %[xx0,zz0]=atdynap(RING0,500,-0.0,0.02);
% 
% [xmaxlist,dplist] = atdynap_om(RING,(0:1:26)*1e-3,1e-6,(-2.4:0.2:2.4)*1e-2,500); 
% %[xmaxlist0,dplist0] = atdynap_om(RING0,(0:1:26)*1e-3,1e-6,(-2.4:0.2:2.4)*1e-2,500); 
% 
% figure
% plot(xx, zz, '.-')
% hold on
% plot(xx0, zz0, 'r.-')
% hold off
% grid on
% set(gca,'fontsize',20);
% xlabel('X [m]');                 % Add labels
% ylabel('Z [m]');
% title('DA (\deltap = 0)')
% addlabel(1, 0, datestr(clock,0))
% %print('Fig1.png','-dpng','-r300')
% 
% figure
% plot(dplist,xmaxlist, '.-')
% hold on
% plot(dplist0,xmaxlist0, 'r.-')
% hold off
% grid on
% set(gca,'fontsize',20);
% title('OFF-momentum DA')
% xlabel('\deltap [%]');                 % Add labels
% ylabel('z [m]');
% %print('Fig2.png','-dpng','-r300')
%%

%atwritem(RING, '/Users/ichaikov/Documents/MATLAB/thomx-mml/machine/THOMX/StorageRing/Lattices/ThomX_016_058_r56_02_chro22')
%atwritem(RING, '/Users/ichaikov/Documents/MATLAB/thomx-mml/machine/THOMX/StorageRing/Lattices/ThomX_017_064_r56_02_chro00')

%%

% atgetfieldvalues(RING,findcells(RING,'FamName','QP31'), 'PolynomB',{1,2})
% atgetfieldvalues(RING,findcells(RING,'FamName','QP4'), 'PolynomB',{1,2})

%%

BareRING = atsetfieldvalues(RING,findcells(RING,'FamName','SX1'), 'PolynomB',{1,3},0);
BareRING = atsetfieldvalues(BareRING,findcells(BareRING,'FamName','SX2'), 'PolynomB',{1,3},0);
BareRING = atsetfieldvalues(BareRING,findcells(BareRING,'FamName','SX3'), 'PolynomB',{1,3},0);

RING = BareRING;
%%

dpp = 0;
[lindata,tune,chrom]=atlinopt(RING,dpp,1:length(RING)+1); 
dispersion=cat(2,lindata.Dispersion)';
beta=cat(1,lindata.beta);
mu=cat(1,lindata.mu);
spos=cat(1,lindata.SPos);

%%

fprintf('#########################" \n')
A=atdata_thomx(RING);
fprintf(' nux    nuz    chix    chiz     bx      bz      Dx      dnux    dnuz    Dx    Emit  \n')
fprintf('%5.2f  %5.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.3f  %6.2f %6.2f  %6.2f  %6.2f  \n',A)
[px ,pz]=get_nudp(0.017, RING);
fprintf('px =  %8.2f  %8.2f   %8.2f   %8.2f   %8.2f  \n',px)
fprintf('pz =  %8.2f  %8.2f   %8.2f   %8.2f   %8.2f  \n',pz)
Q=atdataamp(RING);
fprintf('QA =  %8.2f   %8.2f   %8.2f    \n',Q)
Q=computeRDT(RING, 1, 'tuneshifts');
fprintf('QA =  %8.2f   %8.2f   %8.2f    \n',Q.dnux_dJx,Q.dnux_dJy,Q.dnuy_dJy)
dpp=-0.0;
[xmax]=atdynap_fast(RING,(0:1:40)*1e-3,1e-6,dpp,500);
[zmax]=atdynapz_fast(RING,1e-6,(0:1:20)*1e-3,dpp,500);

QN=get_nuamp([xmax zmax]*0.8,RING,4,dpp);
fprintf('QN =  %8.2f   %8.2f   %8.2f     \n', QN(1)*A(5), QN(2)*A(5), QN(4)*A(6))
fprintf('Xmax DA = %8.2f mm \n',xmax*1e3)
fprintf('Zmax DA = %8.2f mm \n',zmax*1e3)

%figure(100);atplot(RING);
%figure(101);atplot(RING,@plotRDT,'geometric1');
%figure(101);ax=gca;atplot(RING,@plotRDT,'chromatic')

return
%%
DA_beta = load('DA_SXfit1');

[xx,zz]=atdynap(RING,500,0.0,0.02);

RING_quadFF =atsetfieldvalues(RING,find(atgetcells(RING,'Class','Quadrupole')),...
   'PassMethod','QuadMPoleFringePass' );
[xx_QFF,zz_QFF]=atdynap(RING_quadFF,500,-0.0,0.02);


%%
figure(10);
set(gcf,'color','w')
plot(xx*1e3,zz*1e3,'-ob','LineWidth',2,'DisplayName', 'DA in AT FringeBendEntrance/Exit=3 ');
% hold on
% plot(1e3*DA_beta(:,1),1e3*DA_beta(:,2),'r*-','MarkerSize',10,'DisplayName', 'DA in BETA SEXT');
% plot(1e3*xx_QFF,1e3*zz_QFF,'mo-','MarkerSize',7,'DisplayName', 'DA in AT SEXT QUAD FF'); %,'HandleVisibility','off'
% hold off
xlim([-30 30])
ylim([0 15])
grid on
set(gca,'fontsize',20);
xlabel('X [mm]');                 % Add labels
ylabel('Z [mm]');
title('DA (\deltap = 0)')
u = legend('show','Location','NorthEast');
set(u, 'FontSize',16)
addlabel(1, 0, datestr(clock,0))
%print('thomx_DA_FringeBendEntrance3_chro00.png','-dpng','-r300')


%%

dpp = 0;
[lindata,tune,chrom]=atlinopt(RING,dpp,1:length(RING)+1); 
dispersion=cat(2,lindata.Dispersion)';
beta=cat(1,lindata.beta);

%%

rx_bpipe = 20e-3;
rz_bpipe = 14e-3;

dispxinj = dispersion(1,1);
dispzinj = dispersion(1,2);
[dispxmax, ind_x] = max(dispersion(:,1));
[dispzmax, ind_z] = max(dispersion(:,2));

bxinj = beta(1,1);
bzinj = beta(1,2);
% bxdispmax = beta(ind_x,1);
% bzdispmax = beta(ind_z,2);
bxmax = max(beta(:,1));
bzmax = max(beta(:,2));

rx_bpipe_scaled = rx_bpipe / sqrt(bxmax/bxinj)
rz_bpipe_scaled = rz_bpipe / sqrt(bzmax/bzinj)
bpipe_limit = rx_bpipe/dispxmax;



%%

[xmaxlist,dplist] = atdynap_om(RING,(0:2:38)*1e-3,1e-12,(-3:0.2:3)*1e-2,500); 

%%
 
%figure('units','normalized','position',[0.3 0.3 0.45 0.35])
figure(13);
set(gcf,'color','w')
plot(dplist*1e2,xmaxlist*1e3,'-r','LineWidth',2,'DisplayName', 'DA OFF-momentum');
hold on
plot([0 bpipe_limit*1e2], [rx_bpipe_scaled*1e3 0],'k-','DisplayName', 'Scaled vacuum chamber')
plot([0 -bpipe_limit*1e2], [rx_bpipe_scaled*1e3 0],'k-','HandleVisibility','off')
hold off
xlim([-3 3])
ylim([0 max(xmaxlist)]*1.2e3)
grid on
set(gca,'fontsize',20);
%title('OFF-momentum DA')
xlabel('\deltap [%]');                 % Add labels
ylabel('x [mm]');
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
%print('thomx_DAx_offmoment.png','-dpng','-r300')

%%

[zmaxlist,dplist] = atdynapz_om(RING,1e-12,(0:1:20)*1e-3,(-3:0.2:3)*1e-2,500);

figure(11);
set(gcf,'color','w')
plot(dplist*1e2,zmaxlist*1e3,'-b','LineWidth',2,'DisplayName', 'DA OFF-momentum');
hold on
plot([0 max(dplist)*1e2], [rz_bpipe_scaled*1e3 rz_bpipe_scaled*1e3],'k-','DisplayName', 'Scaled vacuum chamber')
plot([0 -max(dplist)*1e2], [rz_bpipe_scaled*1e3 rz_bpipe_scaled*1e3],'k-','HandleVisibility','off')
hold off
xlim([-3 3])
ylim([0 max(zmaxlist)]*1.2e3)
grid on
set(gca,'fontsize',20);
%title('OFF-momentum DA')
xlabel('\deltap [%]');                 % Add labels
ylabel('z [mm]');
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
%print('thomx_DAz_offmoment.png','-dpng','-r300')

%%
rin=[0.003; 0.00; 0.01; 0; 0; 0];
[X,lost]=ringpass(RING,rin,1000); 
figure(50);plot(X(1,:),X(2,:),'.b')
figure(51);plot(X(3,:),X(4,:),'.b')

return