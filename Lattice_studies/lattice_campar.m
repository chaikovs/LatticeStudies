
%%
global GLOBVAL 
GLOBVAL.E0=50E6;
GLOBVAL.LatticeFile='test';

dir ='/Users/ichaikov/Documents/MATLAB/thomx-mml/machine/THOMX/StorageRing/Lattices/';

RING1 = ThomX_017_064_r56_02_chro00();
RING1_multip = ThomX_017_064_r56_02_chro00_multip();
% ring_chro00= ThomX_017_064_r56_02_chro00(); 
RING2= ThomX_016_058_r56_02_chro22(); 
RING2_multip= ThomX_016_058_r56_02_chro22_multip(); 

%%

% RING11=fittune2_alex(RING1, [3.16 1.58], 'QP31', 'QP4');
% RING3=fitchrom_alex(RING11,[2.2 -2.2],'SX2' ,'SX3' );

%RING3 = ThomX_017_056_r56_02_chro22;

%%

atplot(RING1)
atplot(RING2)

atplot(RING1_multip)
atplot(RING2_multip)

%%
fprintf('#########################" \n')
A1=atdata_thomx(RING1);
fprintf(' nux    nuz    chix    chiz     bx      bz      Dx      dnux    dnuz    Dx    Emit  \n')
fprintf('%5.2f  %5.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.3f  %6.2f %6.2f  %6.2f  %6.2f  \n',A1)

A11=atdata_thomx(RING2);
fprintf(' nux    nuz    chix    chiz     bx      bz      Dx      dnux    dnuz    Dx    Emit  \n')
fprintf('%5.2f  %5.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.3f  %6.2f %6.2f  %6.2f  %6.2f  \n',A11)

% A2=atdata_thomx(RING1_multip);
% fprintf(' nux    nuz    chix    chiz     bx      bz      Dx      dnux    dnuz    Dx    Emit  \n')
% fprintf('%5.2f  %5.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.3f  %6.2f %6.2f  %6.2f  %6.2f  \n',A2)
% 
% A2=atdata_thomx(RING2_multip);
% fprintf(' nux    nuz    chix    chiz     bx      bz      Dx      dnux    dnuz    Dx    Emit  \n')
% fprintf('%5.2f  %5.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.3f  %6.2f %6.2f  %6.2f  %6.2f  \n',A2)

%%
[px ,pz]=get_nudp(0.02, RING1);
fprintf('px =  %8.2f   %8.2f   %8.2f   %8.2f  \n',px)
fprintf('pz =  %8.2f   %8.2f   %8.2f   %8.2f  \n',pz)

[px ,pz]=get_nudp(0.02, RING1_multip);
fprintf('px =  %8.2f   %8.2f   %8.2f   %8.2f  \n',px)
fprintf('pz =  %8.2f   %8.2f   %8.2f   %8.2f  \n',pz)

[px ,pz]=get_nudp(0.02, RING2);
fprintf('px =  %8.2f   %8.2f   %8.2f   %8.2f  \n',px)
fprintf('pz =  %8.2f   %8.2f   %8.2f   %8.2f  \n',pz)

[px ,pz]=get_nudp(0.02, RING2_multip);
fprintf('px =  %8.2f   %8.2f   %8.2f   %8.2f  \n',px)
fprintf('pz =  %8.2f   %8.2f   %8.2f   %8.2f  \n',pz)

%%
Q=atdataamp(RING);
fprintf('QA =  %8.2f   %8.2f   %8.2f    \n',Q)

dpp=-0.0;
[xmax]=atdynap_fast(RING,(0:1:40)*1e-3,1e-6,dpp,500);
[zmax]=atdynapz_fast(RING,1e-6,(0:1:20)*1e-3,dpp,500);

QN=get_nuamp([xmax zmax]*0.9,RING,4,dpp);
fprintf('QN =  %8.2f   %8.2f   %8.2f     \n', QN(1)*A(5), QN(2)*A(5), QN(4)*A(6))
fprintf('Xmax DA = %8.2f mm \n',xmax*1e3)
fprintf('Zmax DA = %8.2f mm \n',zmax*1e3)

%%
DA_beta = load('DA_SXfit1');
[xx1,zz1]=atdynap(RING1,500,-0.0,0.02);
[xx1_multip,zz1_multip]=atdynap(RING1_multip,500,-0.0,0.02);
[xx2,zz2]=atdynap(RING2,500,-0.0,0.02);
[xx2_multip,zz2_multip]=atdynap(RING2_multip,500,-0.0,0.02);

%%

indx1=1:length(RING2);    
T=twissring(RING2,0,indx1);
beta=cat(1,T.beta);

%%
rx_bpipe = 20e-3;
rz_bpipe = 14e-3;

bxinj = beta(1,1);
bzinj = beta(1,2);
bxmax = max(beta(:,1));
bzmax = max(beta(:,2));

rx_bpipe_scaled = rx_bpipe / sqrt(bxmax/bxinj)
rz_bpipe_scaled = rz_bpipe / sqrt(bzmax/bzinj)

a=rx_bpipe_scaled; % horizontal radius
b=rz_bpipe_scaled; % vertical radius
x0=0; % x0,y0 ellipse centre coordinates
y0=0;
t=-pi:0.01:pi;
x=x0+a*cos(t);
y=y0+b*sin(t);

a_bpipe=rx_bpipe; % horizontal radius
b_bpipe=rz_bpipe; % vertical radius
x0_bpipe=0; % x0,y0 ellipse centre coordinates
y0_bpipe=0;
t_bpipe=-pi:0.01:pi;
x_bpipe=x0_bpipe+a_bpipe*cos(t_bpipe);
y_bpipe=y0_bpipe+b_bpipe*sin(t_bpipe);

%%
% RING3_quadFF =atsetfieldvalues(RING3,find(atgetcells(RING3,'Class','Quadrupole')),...
%     'PassMethod','QuadMPoleFringePass' );
% [xx3_QFF,zz3_QFF]=atdynap(RING3_quadFF,500,-0.0,0.02);

%%

figure(10);
set(gcf,'color','w')
plot(x*1e3,y*1e3,'k-','DisplayName', 'Scaled vacuum chamber')
hold on
plot(x_bpipe*1e3,y_bpipe*1e3,'k--','DisplayName', 'Vacuum chamber')
plot(xx1*1e3,zz1*1e3,'-ob','LineWidth',2,'DisplayName', 'DA in AT FringeBendEntrance/Exit=3 ');
plot(xx2*1e3,zz2*1e3,'-og','LineWidth',2,'DisplayName', 'DA in AT FringeBendEntrance/Exit=3 ');
plot(xx1_multip*1e3,zz1_multip*1e3,'-om','LineWidth',2,'DisplayName', 'DA in AT FringeBendEntrance/Exit=3 ');
plot(xx2_multip*1e3,zz2_multip*1e3,'-or','LineWidth',2,'DisplayName', 'DA in AT FringeBendEntrance/Exit=3 ');
%plot(1e3*DA_beta(:,1),1e3*DA_beta(:,2),'r*-','MarkerSize',10,'DisplayName', 'DA in BETA SEXT');
%plot(1e3*XX0,1e3*ZZ0,'ko-','MarkerSize',7,'DisplayName', 'DA in AT SEXT'); %,'HandleVisibility','off'
hold off
xlim([-50 50])
ylim([0 30])
grid on
set(gca,'fontsize',20);
xlabel('X [mm]');                 % Add labels
ylabel('Z [mm]');
title('DA (\deltap = 0)')
u = legend('show','Location','NorthEast');
set(u, 'FontSize',16)
addlabel(1, 0, datestr(clock,0))


%%

dpp = 0;
[lindata,tune,chrom]=atlinopt(RING2,dpp,1:length(RING2)+1); 
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

[xmaxlist1,dplist1] = atdynap_om(RING1,(0:2:38)*1e-3,1e-12,(-3:0.2:3)*1e-2,500); 

%%

[xmaxlist2,dplist2] = atdynap_om(RING2,(0:2:38)*1e-3,1e-12,(-3:0.2:3)*1e-2,500); 

%%

[xmaxlist1_multip,dplist1_multip] = atdynap_om(RING1_multip,(0:2:38)*1e-3,1e-12,(-3:0.2:3)*1e-2,500); 

%%

[xmaxlist2_multip,dplist2_multip] = atdynap_om(RING2_multip,(0:2:38)*1e-3,1e-12,(-3:0.2:3)*1e-2,500); 

 %%
%figure('units','normalized','position',[0.3 0.3 0.45 0.35])
figure(13);
set(gcf,'color','w')
plot(dplist1*1e2,xmaxlist1*1e3,'-b','LineWidth',2,'DisplayName', 'DA OFF-momentum 3.17/1.64 Chro 0/0');
hold on
plot(dplist1_multip*1e2,xmaxlist1_multip*1e3,'-r','LineWidth',2,'DisplayName', 'DA OFF-momentum MULTIPOLES 3.17/1.64 Chro 0/0');
plot([0 bpipe_limit*1e2], [rx_bpipe_scaled*1e3 0],'k-','DisplayName', 'Scaled vacuum chamber')
plot([0 -bpipe_limit*1e2], [rx_bpipe_scaled*1e3 0],'k-','HandleVisibility','off')
hold off
xlim([-3 3])
ylim([0 max(xmaxlist1)]*1.4e3)
grid on
set(gca,'fontsize',20);
%title('OFF-momentum DA')
xlabel('\deltap [%]');                 % Add labels
ylabel('x [mm]');
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
%print('OFF-DAhor_Multipoles','-dpng','-r300')

%%

figure(131);
set(gcf,'color','w')
plot(dplist2*1e2,xmaxlist2*1e3,'-b','LineWidth',2,'DisplayName', 'DA OFF-momentum 3.16/1.58 Chro 2/2');
hold on
plot(dplist2_multip*1e2,xmaxlist2_multip*1e3,'-r','LineWidth',2,'DisplayName', 'DA OFF-momentum MULTIPOLES 3.16/1.58 Chro 2/2');
plot([0 bpipe_limit*1e2], [rx_bpipe_scaled*1e3 0],'k-','DisplayName', 'Scaled vacuum chamber')
plot([0 -bpipe_limit*1e2], [rx_bpipe_scaled*1e3 0],'k-','HandleVisibility','off')
hold off
xlim([-3 3])
ylim([0 max(xmaxlist2)]*1.3e3)
grid on
set(gca,'fontsize',20);
%title('OFF-momentum DA')
xlabel('\deltap [%]');                 % Add labels
ylabel('x [mm]');
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
print('OFF-DAhor_Multipoles_WP2','-dpng','-r300')

%%

[zmaxlist1,dplist1] = atdynapz_om(RING1,1e-12,(0:1:20)*1e-3,(-3:0.2:3)*1e-2,500);
[zmaxlist1_multip,dplist1_multip] = atdynapz_om(RING1_multip,1e-12,(0:1:20)*1e-3,(-3:0.2:3)*1e-2,500);

%%
[zmaxlist2,dplist2] = atdynapz_om(RING2,1e-12,(0:1:20)*1e-3,(-3:0.2:3)*1e-2,500);
[zmaxlist2_multip,dplist2_multip] = atdynapz_om(RING2_multip,1e-12,(0:1:20)*1e-3,(-3:0.2:3)*1e-2,500);



%%
figure(14);
set(gcf,'color','w')
plot(dplist1*1e2,zmaxlist1*1e3,'-b','LineWidth',2,'DisplayName', 'DA OFF-momentum 3.17/1.64 Chro 0/0');
hold on
plot(dplist1_multip*1e2,zmaxlist1_multip*1e3,'-r','LineWidth',2,'DisplayName', 'DA OFF-momentum MULTIPOLES 3.17/1.64 Chro 0/0');
plot([0 max(dplist1)*1e2], [rz_bpipe_scaled*1e3 rz_bpipe_scaled*1e3],'k-','DisplayName', 'Scaled vacuum chamber')
plot([0 -max(dplist1)*1e2], [rz_bpipe_scaled*1e3 rz_bpipe_scaled*1e3],'k-','HandleVisibility','off')
hold off
xlim([-3 3])
ylim([0 max(zmaxlist1)]*1.4e3)
grid on
set(gca,'fontsize',20);
%title('OFF-momentum DA')
xlabel('\deltap [%]');                 % Add labels
ylabel('z [mm]');
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
%print('OFF-DAvert_Multipoles','-dpng','-r300')

%%

figure(141);
set(gcf,'color','w')
plot(dplist2*1e2,zmaxlist2*1e3,'-b','LineWidth',2,'DisplayName', 'DA OFF-momentum 3.16/1.58 Chro 2/2');
hold on
plot(dplist2_multip*1e2,zmaxlist2_multip*1e3,'-r','LineWidth',2,'DisplayName', 'DA OFF-momentum MULTIPOLES 3.16/1.58 Chro 2/2');
plot([0 max(dplist2)*1e2], [rz_bpipe_scaled*1e3 rz_bpipe_scaled*1e3],'k-','DisplayName', 'Scaled vacuum chamber')
plot([0 -max(dplist2)*1e2], [rz_bpipe_scaled*1e3 rz_bpipe_scaled*1e3],'k-','HandleVisibility','off')
hold off
xlim([-3 3])
ylim([0 max(zmaxlist2)]*1.3e3)
grid on
set(gca,'fontsize',20);
%title('OFF-momentum DA')
xlabel('\deltap [%]');                 % Add labels
ylabel('z [mm]');
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
print('OFF-DAvert_Multipoles_WP22','-dpng','-r300')




