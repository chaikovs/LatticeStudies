
clc; close all; clear all;


%% Initial lattices

thomx_ring=ThomX_017_064_r56_03_chro11;



%% As you wrote in the email, I am getting NaN if x=y=0.3 which indicates perhpas that particles are lost ???

% Check out this :
% rin=[0.03; 0.00; 0.03; 0; 0; 0];
% [X,lost]=ringpass(RING,rin,1000,'reuse'); 
% figure(50)
% plot(X(1,:),X(2,:),'.b')


rin0=[0.001; 0.00; 0.001; 0; 0; 0];
rin=[0.005; 0.00; 0.005; 0; 0; 0];
rin2=[0.01; 0.00; 0.01; 0; 0; 0];
rin3=[0.013; 0.00; 0.013; 0; 0; 0];
rin4=[0.015; 0.00; 0.015; 0; 0; 0];

% rin0p=[0.005; 0.00; 0.005; 0; -0.01; 0];
% rinp=[0.03; 0.00; 0.03; 0; 0; 0];
% rin2p=[0.01; 0.00; 0.01; 0; 0.01; 0.0];
% rin3p=[0.013; 0.00; 0.013; 0; 0.015; 0];
% rin4p=[0.015; 0.00; 0.015; 0; 0.018; 0];

rin0p=[0.001; 0.00; 0.001; 0; -0.01; 0];
rinp=[0.001; 0.00; 0.001; 0; 0; 0];
rin2p=[0.001; 0.00; 0.001; 0; 0.01; 0.0];
rin3p=[0.001; 0.00; 0.001; 0; 0.015; 0];
rin4p=[0.001; 0.00; 0.001; 0; 0.018; 0];

[X0,lost0]=ringpass(thomx_ring,rin0,1000);
[X,lost]=ringpass(thomx_ring,rin,1000); 
[X2,lost2]=ringpass(thomx_ring,rin2,1000);
[X3,lost3]=ringpass(thomx_ring,rin3,1000);
[X4,lost4]=ringpass(thomx_ring,rin4,1000);

[X0p,lost0p]=ringpass(thomx_ring,rin0p,1000);
[Xp,lostp]=ringpass(thomx_ring,rinp,1000); 
[X2p,lost2p]=ringpass(thomx_ring,rin2p,1000);
[X3p,lost3p]=ringpass(thomx_ring,rin3p,1000);
[X4p,lost4p]=ringpass(thomx_ring,rin4p,1000);


figure(51)
set(gca,'fontsize',16);
plot(X2(1,:),X2(2,:),'.m','DisplayName', 'x=y=0.01 \deltaP=0')
hold on
plot(X(1,:),X(2,:),'.g','DisplayName', 'x=y=0.005 \deltaP=0')
plot(X0(1,:),X0(2,:),'.b','DisplayName', 'x=y=0.001 \deltaP=0')
plot(X3(1,:),X3(2,:),'*k','DisplayName', 'x=y=0.013 \deltaP=0')
plot(X4(1,:),X4(2,:),'or','DisplayName', 'x=y=0.015 \deltaP=0')
%xlim([-3 3])
%ylim([0 max(xmaxlist)]*1.2e3)
grid on
xlabel('x [m]');                 % Add labels
ylabel('xp [rad]');
u = legend('show','Location','SouthEast');
set(u,'FontSize',14)
%print('thomx_4Dtracking_r56_0p4.png','-dpng','-r300')

figure(52)
plot(X2p(5,:),X2p(6,:),'.m','DisplayName', '\deltaP=0.01 x=y=0.001')
hold on
plot(Xp(5,:),Xp(6,:),'dg','DisplayName', '\deltaP=0.0 x=y=0.001')
plot(X0p(5,:),X0p(6,:),'.b','DisplayName', '\deltaP=-0.01 x=y=0.001')
plot(X3p(5,:),X3p(6,:),'*k','DisplayName', '\deltaP=0.015 x=y=0.001')
plot(X4p(5,:),X4p(6,:),'or','DisplayName', '\deltaP=0.018 x=y=0.001')
%xlim([-3 3])
%ylim([0 max(xmaxlist)]*1.2e3)
grid on
set(gca,'fontsize',20);
xlabel('\deltap');                 % Add labels
ylabel('c\tau [m]');
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
%print('thomx_6Dtracking_r56_0p4.png','-dpng','-r300')


%% DA off-momentum


thomx_ring2=atcavityoff(thomx_ring);

%%
[XX,ZZ]   = atdynap(thomx_ring, 500,0,0.02); 
[XX2,ZZ2]   = atdynap(thomx_ring2, 500,0,0.02); 

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
plot(XX,ZZ,'ro-','MarkerSize',7,'LineWidth',2,'DisplayName', 'RF ON'); %,'HandleVisibility','off'
hold on; 
plot(XX2,ZZ2,'bo-','MarkerSize',7,'LineWidth',2,'DisplayName', 'RF OFF'); %,'HandleVisibility','off'
xlabel('x [m]')
ylabel('z [m]')
set(gca,'FontSize',20)
set(gcf,'color','w')
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
xlim([-0.04 0.04])
ylim([0 0.023])
addlabel(1, 0, datestr(clock,0))
print('DAx_onmoment_thomx_r56_0p3.png','-dpng','-r300')


%%

[xmaxlist,dplist] = atdynap_om(thomx_ring,(0:2:38)*1e-3,1e-12,(-3:0.2:3)*1e-2,500); 


%%
[xmaxlist2,dplist2] = atdynap_om(thomx_ring2,(0:2:38)*1e-3,1e-12,(-3:0.2:3)*1e-2,500); 


%%
 
%figure('units','normalized','position',[0.3 0.3 0.45 0.35])
figure(13);
set(gcf,'color','w')
plot(dplist*1e2,xmaxlist*1e3,'-r','LineWidth',2,'DisplayName', 'DA OFF-momentum with RF');
hold on
plot(dplist2*1e2,xmaxlist2*1e3,'-b','LineWidth',2,'DisplayName', 'DA OFF-momentum without RF');
%plot([0 bpipe_limit*1e2], [rx_bpipe_scaled*1e3 0],'k-','DisplayName', 'Scaled vacuum chamber')
%plot([0 -bpipe_limit*1e2], [rx_bpipe_scaled*1e3 0],'k-','HandleVisibility','off')
hold off
%xlim([-3 3])
%ylim([0 max(xmaxlist)]*1.2e3)
grid on
set(gca,'fontsize',20);
%title('OFF-momentum DA')
xlabel('\deltap [%]');                 % Add labels
ylabel('x [mm]');
u = legend('show','Location','SouthWest');
set(u,'FontSize',14)
print('DAx_offmoment_thomx_r56_0p3.png','-dpng','-r300')



