

%%
% 
% DeltaRF = [ -0.2 -0.1 -0.05 -0.02 0.00001 0.02 0.05 0.1 0.2]*1e6;
% 
% %DeltaRF = getrf('Physics') * getmcf * dpdelta 
% RF3 = getrf('Physics');
% RF0=RF3(1);
% %dp = DeltaRF / (RF0*MCF)
% dp = DeltaRF ./ (RF0*getmcf())
% 
% for irf = 1:length(DeltaRF)
%     
% [Chromaticity,Tune] = modelchro(DeltaRF(irf))
% 
% chrx(irf) = Chromaticity(1,1);
% chry(irf) = Chromaticity(2,1);
% 
% nux(irf) = Tune(1,1);
% nuy(irf) = Tune(2,1);
% 
% end
% 
% plot(dp,nuy,'.-' )
%%
%RING = ThomX_017_064_r56_02_chro00;

 RING = thomx_4NEWpass;
% RING=atsetfieldvalues(RING,findcells(RING,'PassMethod','BndMPoleSymplecticNew4Pass'),...
%     'PassMethod','BndMPoleSymplectic4Pass');
% 
% ind = findcells(RING,'PassMethod','BndMPoleSymplectic4Pass')
% 
% for i = 1:8
% RING{ind(i)}.FringeBendEntrance = 3;
% RING{ind(i)}.FringeBendExit = 3;
% end

tuneDPx = [];
tuneDPz = [];
tuneDPx_NEW = [];
tuneDPz_NEW = [];
tuneDPx_NEW_wo_FF = [];
tuneDPz_NEW_wo_FF = [];
dp = -0.05:0.002:0.05;
for idp = dp
[TD,t,c] = atlinopt(RING,idp,1);
tuneDPx = [tuneDPx t(1)];
tuneDPz = [tuneDPz t(2)];

end
[TD,t0,c] = atlinopt(RING,0,1);

figure
set(gca,'FontSize',18)
plot(dp,tuneDPx - t0(1) ,'r.-','MarkerSize',9 )
hold on
plot(dp,tuneDPz - t0(2) ,'b.-','MarkerSize',9 )
hold off


%%
ring_NEW = thomx_4NEWpass();
ring = thomx_4pass();
ring_NEW_wo_FF = thomx_4NEW_wo_FF();

ring_linear=atsetfieldvalues(ring,findcells(ring,'PassMethod','BndMPoleSymplecticNew4Pass'),...
    'PassMethod','BendLinearPass');

tuneDPx = [];
tuneDPz = [];
tuneDPx_NEW = [];
tuneDPz_NEW = [];
tuneDPx_NEW_wo_FF = [];
tuneDPz_NEW_wo_FF = [];

dp = -0.05:0.005:0.05;
for idp = dp
[TD,t,c] = atlinopt(ring,idp,1);
tuneDPx = [tuneDPx t(1)];
tuneDPz = [tuneDPz t(2)];

end

dp_NEW = -0.06:0.005:0.06;
for idp = dp_NEW
[TD_NEW,t_NEW,c_NEW] = atlinopt(ring_NEW,idp,1);
tuneDPx_NEW = [tuneDPx_NEW t_NEW(1)];
tuneDPz_NEW = [tuneDPz_NEW t_NEW(2)];

end

dp_NEW_wo_FF = -0.05:0.005:0.05;
for idp = dp_NEW_wo_FF
[TD_NEW_wo_FF,t_NEW_wo_FF,c_NEW_wo_FF] = atlinopt(ring_NEW_wo_FF,idp,1);
tuneDPx_NEW_wo_FF = [tuneDPx_NEW_wo_FF t_NEW_wo_FF(1)];
tuneDPz_NEW_wo_FF = [tuneDPz_NEW_wo_FF t_NEW_wo_FF(2)];

end

tuneDPx_linear = [];
tuneDPz_linear = [];
dp_linear = -0.06:0.005:0.06;
for idp = dp_linear
[TD_linear,t_linear,c_linear] = atlinopt(ring_linear,idp,1);
tuneDPx_linear = [tuneDPx_linear t_NEW(1)];
tuneDPz_linear = [tuneDPz_linear t_NEW(2)];

end

[TD,t0,c] = atlinopt(ring_NEW,0,1);
%% BETA anal dypole

data_analdip_nux = load('nux_analdipole');
data_analdip_nuz = load('nuz_analdipole');
data_dip_nux = load('nux_ldipole');
data_dip_nuz = load('nuz_ldipole');


% figure
% subplot 211
% set(gca,'FontSize',18)
% plot(data_analdip_nux(:,1),data_analdip_nux(:,2) - t0(1) ,'rd-' )
% hold on
% plot(data_dip_nux(:,1),data_dip_nux(:,2) - t0(1) ,'r*-' )
% ylabel('Horizontal tune difference')
% title(['Tune shift with dp/p ' ' Working Point ' num2str(round(t0,4))])
% subplot 212
% set(gca,'FontSize',18)
% plot(data_analdip_nuz(:,1),data_analdip_nuz(:,2) - t0(2) ,'bd-' )
% hold on
% plot(data_dip_nuz(:,1),data_dip_nuz(:,2) - t0(2) ,'b*-' )
% hold off
% xlabel(' Momentum diviation dp/p')
% ylabel('Vertical tune difference')
%% Alex code

data_alex_FF = load('code-compar/CODE_TDR_0.17_0.64_r56_0.2_sx_Dff412_DipMagnL_chro00-NUDP.mat');
data_alex = load('code-compar/CODE_TDR_0.17_0.64_r56_0.2_sx_412_DipMagnL_chro00-NUDP.mat');



%% With fringe


figure
subplot 211
set(gca,'FontSize',18)
plot(dp_NEW,tuneDPx_NEW - t0(1) ,'r.-','MarkerSize',9 ,'DisplayName', 'AT (BndMPoleSymplecticNew4Pass)')
hold on
plot(dp,tuneDPx - t0(1) ,'rs-','MarkerSize',4 ,'DisplayName', 'AT (BndMPoleSymplectic4Pass)')
%plot(data_analdip_nux(:,1),data_analdip_nux(:,2) - t0(1) ,'kd-','MarkerSize',4,'DisplayName', 'BETA (analytical dipole)')
plot(data_dip_nux(:,1),data_dip_nux(:,2) - t0(1) ,'k*-','MarkerSize',4,'DisplayName', 'BETA (fringing field)' )
plot(data_alex_FF.dp2,data_alex_FF.nux2- t0(1)  ,'m*-','MarkerSize',4,'DisplayName', 'Alex code (fringing field)' )
ylabel('Horizontal tune difference')
title(['Tune shift with dp/p. With fringing field. ' ' Working Point: ' num2str(round(t0,4))])
legend('show','Location','SouthEast');
subplot 212
set(gca,'FontSize',18)
plot(dp_NEW,tuneDPz_NEW -t0(2) ,'b.-','MarkerSize',9,'DisplayName', 'AT (BndMPoleSymplecticNew4Pass)' )
hold on
plot(dp,tuneDPz -t0(2) ,'bs-','MarkerSize',4,'DisplayName', 'AT (BndMPoleSymplectic4Pass)' )
%plot(data_analdip_nuz(:,1),data_analdip_nuz(:,2) - t0(2) ,'kd-','MarkerSize',4,'DisplayName', 'BETA (analytical dipole)' )
plot(data_dip_nuz(:,1),data_dip_nuz(:,2) - t0(2) ,'k*-','MarkerSize',4 ,'DisplayName', 'BETA (fringing field)')
plot(data_alex_FF.dp2,data_alex_FF.nuz2 - t0(2)  ,'m*-','MarkerSize',4,'DisplayName', 'Alex code (fringing field)' )
hold off
ylim([-0.4 0.4])
xlabel(' Momentum deviation dp/p')
ylabel('Vertical tune difference')
legend('show','Location','SouthEast');
print('thomx_tunediff_dpp_FF.png','-dpng','-r300')


%% Without fringe


figure
subplot 211
set(gca,'FontSize',18)
plot(dp_NEW,tuneDPx_NEW - t0(1) ,'r.-','MarkerSize',9 ,'DisplayName', 'AT with FF (BndMPoleSymplecticNew4Pass)')
hold on
plot(dp_NEW_wo_FF,tuneDPx_NEW_wo_FF - t0(1) ,'g.-','MarkerSize',9 ,'DisplayName', 'AT wo FF (BndMPoleSymplecticNew4Pass)')
%plot(dp,tuneDPx - t0(1) ,'rs-','MarkerSize',4 ,'DisplayName', 'AT (BndMPoleSymplectic4Pass)')
plot(data_analdip_nux(:,1),data_analdip_nux(:,2) - t0(1) ,'kd-','MarkerSize',4,'DisplayName', 'BETA (analytical dipole)')
%plot(data_dip_nux(:,1),data_dip_nux(:,2) - t0(1) ,'k*-','MarkerSize',4,'DisplayName', 'BETA (fringing field)' )
plot(data_alex.dp2,data_alex.nux2- t0(1)  ,'m*-','MarkerSize',4,'DisplayName', 'Alex code' )
ylabel('Horizontal tune difference')
title(['Tune shift with dp/p. Without fringing field. ' ' Working Point: ' num2str(round(t0,4))])
legend('show','Location','SouthEast');

subplot 212
set(gca,'FontSize',18)
plot(dp_NEW,tuneDPz_NEW -t0(2) ,'b.-','MarkerSize',9,'DisplayName', 'AT with FF (BndMPoleSymplecticNew4Pass)' )
hold on
plot(dp_NEW_wo_FF,tuneDPz_NEW_wo_FF - t0(2) ,'g.-','MarkerSize',9 ,'DisplayName', 'AT wo FF (BndMPoleSymplecticNew4Pass)')
%plot(dp,tuneDPz -t0(2) ,'bs-','MarkerSize',4,'DisplayName', 'AT (BndMPoleSymplectic4Pass)' )
plot(data_analdip_nuz(:,1),data_analdip_nuz(:,2) - t0(2) ,'kd-','MarkerSize',4,'DisplayName', 'BETA (analytical dipole)' )
%plot(data_dip_nuz(:,1),data_dip_nuz(:,2) - t0(2) ,'k*-','MarkerSize',4 ,'DisplayName', 'BETA (fringing field)')
plot(data_alex.dp2,data_alex.nuz2 - t0(2)  ,'m*-','MarkerSize',4,'DisplayName', 'Alex code' )
hold off
ylim([-0.4 0.6])
xlabel(' Momentum deviation dp/p')
ylabel('Vertical tune difference')
legend('show','Location','NorthWest');
print('thomx_tunediff_dpp.png','-dpng','-r300')


%% plot for the tune by itself

figure
subplot 211
set(gca,'FontSize',18)
plot(dp_NEW,tuneDPx_NEW  ,'r.-','MarkerSize',9 ,'DisplayName', 'AT (BndMPoleSymplecticNew4Pass)')
hold on
plot(dp,tuneDPx  ,'rs-','MarkerSize',4 ,'DisplayName', 'AT (BndMPoleSymplectic4Pass)')
%plot(dp_linear,tuneDPx_linear  ,'ms-','MarkerSize',4 ,'DisplayName', 'AT (BendLinearPass)')
plot(data_analdip_nux(:,1),data_analdip_nux(:,2) ,'kd-','MarkerSize',4,'DisplayName', 'BETA (analytical dipole)')
plot(data_dip_nux(:,1),data_dip_nux(:,2)  ,'k*-','MarkerSize',4,'DisplayName', 'BETA (fringing field)' )
plot(data_alex.dp2,data_alex.nux2  ,'g*-','MarkerSize',4,'DisplayName', 'Alex (fringing field)' )
ylabel('Horizontal tune ')
title(['Tune shift with dp/p. ' ' Working Point: ' num2str(round(t0,4))])
legend('show','Location','SouthEast');
subplot 212
set(gca,'FontSize',18)
plot(dp_NEW,tuneDPz_NEW  ,'b.-','MarkerSize',9,'DisplayName', 'AT (BndMPoleSymplecticNew4Pass)' )
hold on
plot(dp,tuneDPz ,'bs-','MarkerSize',4,'DisplayName', 'AT (BndMPoleSymplectic4Pass)' )
%plot(dp_linear,tuneDPz_linear  ,'ms-','MarkerSize',4 ,'DisplayName', 'AT (BendLinearPass)')
plot(data_analdip_nuz(:,1),data_analdip_nuz(:,2)  ,'kd-','MarkerSize',4,'DisplayName', 'BETA (analytical dipole)' )
plot(data_dip_nuz(:,1),data_dip_nuz(:,2) ,'k*-','MarkerSize',4 ,'DisplayName', 'BETA (fringing field)')
plot(data_alex.dp2,data_alex.nuz2  ,'g*-','MarkerSize',4,'DisplayName', 'Alex (fringing field)' )
hold off
ylim([0.1 1.1])
xlabel(' Momentum deviation dp/p')
ylabel('Vertical tune ')
legend('show','Location','SouthEast');
print('thomx_tune_dpp.png','-dpng','-r300')


%%

r_QP1 = findcells(ring,'FamName','QP1')
RING_QP1 = findcells(RING,'FamName','QP1')

% indQCor=find(atgetcells(RING,'Class','Quadrupole'))'

atgetfieldvalues(ring,r_QP1,'PolynomB',{1,2})
atgetfieldvalues(RING,RING_QP1,'PolynomB',{1,2})

