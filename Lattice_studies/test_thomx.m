clear
global GLOBVAL 
GLOBVAL.E0=50E6;
GLOBVAL.LatticeFile='test';

% dir ='/home/sources/physmach/loulergue/work/matlab/Simu;
dir ='/Users/ichaikov/Documents/MATLAB/thomx-mml/machine/THOMX/StorageRing/Lattices/BETA/';
%input_file =[dir 'cls/ThomX-30-V5/Good/TDR_0.17_0.64_r56_0.2_sx_Dff41.2_DipMagnL_chro_1_1.str'];
%input_file =[dir 'cls/ThomX-30-V5/Good-old/TDR_0.17_0.64_r56_0.2_sx_Dff.str'];
%input_file =[dir 'cls/ThomX-30-V5/Iryna2/TDR_0.17_0.64_r56_0.2_sx_Dff41.2_FF_chro00'];
input_file =[dir 'TDR_0.17_0.64_r56_0.2_sx_Dff41.2_FF_chro00_SXfit1'];
%input_file =[dir 'cls/ThomX-30-V5/divers/TDR_0.17_0.64_r56_0.2_sx_Dff_tracy.str'];
RING=atreadbeta_thomx_alex(input_file);

%RING=scalesext(RING,'SX1',0.1);
% RING=scalesext(RING,'SX2',1);
% RING=scalesext(RING,'SX3',1.3);
%RING=fitchrom_alex(RING,[0. -0.0],'SX2' ,'SX3' );

% dnu=[0.0 0.2];
% for k=1:5
% RING=fittune2_alex(RING, [(3.17+k*dnu(1)/5) (1.64+k*dnu(2)/5)], 'QP3', 'QP41'); 
% end
 atplot(RING)
[l,t,c] = atlinopt(RING,0,1);
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
fprintf('#########################" \n')
A=atdata_thomx(RING);
fprintf(' nux    nuz    chix    chiz     bx      bz      Dx      dnux    dnuz    Dx    Emit  \n')
fprintf('%5.2f  %5.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.3f  %6.2f %6.2f  %6.2f  %6.2f  \n',A)
[px ,pz]=get_nudp(0.05, RING);
fprintf('px =  %8.2f   %8.2f   %8.2f   %8.2f  \n',px)
fprintf('pz =  %8.2f   %8.2f   %8.2f   %8.2f  \n',pz)
Q=atdataamp(RING);
fprintf('QA =  %8.2f   %8.2f   %8.2f    \n',Q)
Q=computeRDT(RING, 1, 'tuneshifts');
fprintf('QA =  %8.2f   %8.2f   %8.2f    \n',Q.dnux_dJx,Q.dnux_dJy,Q.dnuy_dJy)
dpp=-0.0;
[xmax]=atdynap_fast(RING,(0:1:40)*1e-3,1e-6,dpp,500);
[zmax]=atdynapz_fast(RING,1e-6,(0:1:20)*1e-3,dpp,500);
QN=get_nuamp([xmax zmax]*0.9,RING,4,dpp);
fprintf('QN =  %8.2f   %8.2f   %8.2f     \n', QN(1)*A(5), QN(2)*A(5), QN(4)*A(6))
fprintf('Xmax DA = %8.2f mm \n',xmax*1e3)
fprintf('Zmax DA = %8.2f mm \n',zmax*1e3)

figure(100);atplot(RING);
%figure(101);atplot(RING,@plotRDT,'geometric1');
%figure(101);ax=gca;atplot(RING,@plotRDT,'chromatic')

return
%%
DA_beta = load('DA_SXfit1');

[xx,zz]=atdynap(RING,500,-0.0,0.02);

%%
figure(10);
set(gcf,'color','w')
plot(xx*1e3,zz*1e3,'-ob','LineWidth',2,'DisplayName', 'DA in AT FringeBendEntrance/Exit=3 ');
hold on
plot(1e3*DA_beta(:,1),1e3*DA_beta(:,2),'r*-','MarkerSize',10,'DisplayName', 'DA in BETA SEXT');
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
%print('thomx_DA_FringeBendEntrance3_chro00.png','-dpng','-r300')

%%
[xmaxlist,dplist] = atdynap_om(RING,(0:2:40)*1e-3,1e-12,(-3:0.2:3)*1e-2,500); 
figure(13);
set(gcf,'color','w')
set(gca,'fontsize',16');
plot(dplist*1e2,xmaxlist*1e3,'-b','LineWidth',2);
xlim([-3 3])
ylim([0 max(xmaxlist)]*1.2e3)
grid on
xlabel('dp [%]');                 % Add labels
ylabel('X [mm]');

%%
[zmaxlist,dplist] = atdynapz_om(RING,1e-12,(0:1:40)*1e-3,(-3:0.2:3)*1e-2,500);
figure(11);
set(gcf,'color','w')
set(gca,'fontsize',16');
plot(dplist*1e2,zmaxlist*1e3,'-b','LineWidth',2);
xlim([-3 3])
ylim([0 max(zmaxlist)]*1.2e3)
grid on
xlabel('dp [%]');                 % Add labels
ylabel('Z [mm]');

%%
rin=[0.003; 0.00; 0.01; 0; 0; 0];
[X,lost]=ringpass(RING,rin,100); 
figure(50);plot(X(1,:),X(2,:),'.b')
figure(51);plot(X(3,:),X(4,:),'.b')

return