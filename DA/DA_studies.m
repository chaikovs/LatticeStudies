
%%
% global THERING
% indx=1:length(THERING);  
% T=twissring(THERING,0,indx);   
% beta=cat(1,T.beta);
%[bx by] = modelbeta;
%%

ring0= ThomX_017_064_r56_02(); 
ring = ThomX_017_064_r56_02_chro00();
ring_multip = ThomX_017_064_r56_02_chro00_multip_AT2();
% ring_chro00= ThomX_017_064_r56_02_chro00(); 
%ring2= ThomX_016_058_r56_02_chro22(); 
%ring2_multip= ThomX_016_058_r56_02_chro22_multip(); 

ring_WOquadFF =atsetfieldvalues(ring,find(atgetcells(ring,'Class','Quadrupole')),...
    'PassMethod','StrMPoleSymplectic4Pass' );


%%
indx=1:length(ring);    
T=twissring(ring,0,indx);
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

DA_beta = load('DA_SXfit1');
[XX0,ZZ0]   = atdynap(ring0, 500,0,0.02);
[XX,ZZ]   = atdynap(ring, 500,0,0.02); 
[XX_multip,ZZ_multip]   = atdynap(ring_multip, 500,0,0.02); 

[XX_WOquadFF,ZZ_WOquadFF]   = atdynap(ring_WOquadFF, 500,0,0.02); 

% [XX2,ZZ2]   = atdynap(ring2, 500,0,0.02); 
% [XX2_multip,ZZ2_multip]   = atdynap(ring2_multip, 500,0,0.02);

[l,t,c] = atlinopt(ring,0,1);
[l_mp,t_mp,c_mp] = atlinopt(ring_multip,0,1);

%%

ring_quadFF =atsetfieldvalues(ring,find(atgetcells(ring,'Class','Quadrupole')),...
    'PassMethod','QuadMPoleFringePass' );

ring_multip_quadFF =atsetfieldvalues(ring_multip,find(atgetcells(ring_multip,'Class','Quadrupole')),...
    'PassMethod','QuadMPoleFringePass' );

% ring2_quadFF =atsetfieldvalues(ring2,find(atgetcells(ring2,'Class','Quadrupole')),...
%     'PassMethod','QuadMPoleFringePass' );
% 
% ring2_multip_quadFF =atsetfieldvalues(ring2_multip,find(atgetcells(ring2_multip,'Class','Quadrupole')),...
%     'PassMethod','QuadMPoleFringePass' );

[XX_quadFF,ZZ_quadFF]   = atdynap(ring_quadFF, 500,0.0,0.02); 
[XX_multip_quadFF,ZZ_multip_quadFF]   = atdynap(ring_multip_quadFF, 500,0.0,0.02); 

% [XX2_quadFF,ZZ2_quadFF]   = atdynap(ring2_quadFF, 500,0.0,0.02); 
% [XX2_multip_quadFF,ZZ2_multip_quadFF]   = atdynap(ring2_multip_quadFF, 500,0.0,0.02); 

DA_beta_quadFF = load('DA_SXfit1_Qff');

%%

% ring_temp=scalesext(ring,'SX1',0.1);
% %ring0=atfitchrom(ring_temp,[0.0 -0.0],'SX2','SX3');
% ring0=fitchrom_alex(ring_temp,[0. -0.0],'SX2' ,'SX3' );
% [XX0,ZZ0]  = atdynap(ring0, 500,0,0.02); 
% 
% %[TD, tunes, chrom] = twissring(ring, 0, 1:length(ring)+1,'chrom', 1e-8);
% 
% [l0,t0,c0] = atlinopt(ring0,0,1);


%%

% S_sx1 = atgetfieldvalues(ring,findcells(ring,'FamName','SX1'), 'PolynomB',{1,3});
% S_sx2 = atgetfieldvalues(ring,findcells(ring,'FamName','SX2'), 'PolynomB',{1,3});
% S_sx3 = atgetfieldvalues(ring,findcells(ring,'FamName','SX3'), 'PolynomB',{1,3});
% 
% 
% S0_sx1 = atgetfieldvalues(ring0,findcells(ring0,'FamName','SX1'), 'PolynomB',{1,3});
% S0_sx2 = atgetfieldvalues(ring0,findcells(ring0,'FamName','SX2'), 'PolynomB',{1,3});
% S0_sx3 = atgetfieldvalues(ring0,findcells(ring0,'FamName','SX3'), 'PolynomB',{1,3});
% 
% fprintf('SX1 = %8.2f ; SX10 = %8.2f  \n', S_sx1(1), S0_sx1(1))
% fprintf('SX2 = %8.2f ; SX20 = %8.2f  \n', S_sx2(1), S0_sx2(1))
% fprintf('SX3 = %8.2f ; SX30 = %8.2f  \n', S_sx3(1), S0_sx3(1))

%%

% figure('units','normalized','position',[0.3 0.3 0.45 0.35])
% plot(x,y,'k-','DisplayName', 'Scaled vacuum chamber')
% hold on; 
% plot(x_bpipe,y_bpipe,'k--','DisplayName', 'Vacuum chamber')
% plot(XX,ZZ,'bo-','MarkerSize',7,'LineWidth',2,'DisplayName', ['DA in AT SEXT Chro '  num2str(c)]); %,'HandleVisibility','off'
% plot(DA_beta(:,1),DA_beta(:,2),'r*-','MarkerSize',10,'LineWidth',2,'DisplayName', 'DA in BETA SEXT Chro 00');
% %plot(XX0,ZZ0,'mo-','MarkerSize',7,'DisplayName', 'DA in AT SEXT Chro 00'); %,'HandleVisibility','off'
% plot(XX_chro00,ZZ_chro00,'mo-','MarkerSize',7,'LineWidth',2,'DisplayName', ['DA in AT SEXT Chro ' num2str(c_chro00)]); %,'HandleVisibility','off'
% xlabel('x [m]')
% ylabel('z [m]')
% set(gca,'FontSize',20)
% set(gcf,'color','w')
% u = legend('show','Location','NorthEast');
% set(u,'FontSize',14)
% xlim([-0.05 0.05])
% ylim([0 0.025])
% addlabel(1, 0, datestr(clock,0))
% print('DA_comparizon_BETA_AT_chro00','-dpng','-r300')

%%

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
plot(x,y,'k-','DisplayName', 'Scaled vacuum chamber')
hold on; 
plot(x_bpipe,y_bpipe,'k--','DisplayName', 'Vacuum chamber')
%plot(XX,ZZ,'bo-','MarkerSize',7,'LineWidth',2,'DisplayName', ['DA in AT SEXT Chro '  num2str(c)]); %,'HandleVisibility','off'
plot(DA_beta(:,1),DA_beta(:,2),'r*-','MarkerSize',10,'LineWidth',2,'DisplayName', 'DA in BETA: SEXT Chro 00');
%plot(XX0,ZZ0,'mo-','MarkerSize',7,'DisplayName', 'DA in AT SEXT Chro 00'); %,'HandleVisibility','off'
plot(XX_WOquadFF,ZZ_WOquadFF,'mo-','MarkerSize',7,'LineWidth',2,'DisplayName', ['DA in AT: SEXT Chro (fitted) ' num2str(c)]); %,'HandleVisibility','off'
xlabel('x [m]')
ylabel('z [m]')
set(gca,'FontSize',20)
set(gcf,'color','w')
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
xlim([-0.05 0.05])
ylim([0 0.025])
addlabel(1, 0, datestr(clock,0))
print('DA_BETA_AT_chro00','-dpng','-r300')


%%


figure('units','normalized','position',[0.3 0.3 0.45 0.35])
plot(x,y,'k-','DisplayName', 'Scaled vacuum chamber')
hold on; 
plot(x_bpipe,y_bpipe,'k--','DisplayName', 'Vacuum chamber')
%plot(XX,ZZ,'bo-','MarkerSize',7,'DisplayName', 'DA in AT SEXT'); %,'HandleVisibility','off'
%plot(DA_beta(:,1),DA_beta(:,2),'ro-','MarkerSize',7,'DisplayName', 'DA in BETA SEXT');
plot(XX_WOquadFF,ZZ_WOquadFF,'mo-','MarkerSize',7,'LineWidth',2,'DisplayName', ['DA in AT: SEXT Chro (fitted) ' num2str(c)]); %,'HandleVisibility','off'
plot(DA_beta_quadFF(:,1),DA_beta_quadFF(:,2),'r*-','MarkerSize',10,'LineWidth',2,'DisplayName', 'DA in BETA: SEXT + QUAD Fringing Field');
plot(XX,ZZ,'b*-','MarkerSize',10,'LineWidth',2,'DisplayName', 'DA in AT: SEXT + QUAD Fringing Field');
%plot(XX0,ZZ0,'b*-','DisplayName', 'DA in AT Chro 0');
% plot(XXq,ZZq,'mo-');
% plot(XXq0,ZZq0,'m*-');
xlabel('x [m]')
ylabel('z [m]')
set(gca,'FontSize',20)
set(gcf,'color','w')
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
xlim([-0.05 0.05])
ylim([0 0.025])
print('DA_QuadFF','-dpng','-r300')

%%

DA_beta_multip = load('DA_SXfit1_DipQuadSextM');
DA_beta_multip_Qff = load('DA_SXfit1_Qff_DipQuadSextM');

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
plot(x,y,'k-','DisplayName', 'Scaled vacuum chamber')
hold on; 
plot(x_bpipe,y_bpipe,'k--','DisplayName', 'Vacuum chamber')
plot(XX_WOquadFF,ZZ_WOquadFF,'mo-','MarkerSize',7,'LineWidth',2,'DisplayName', ['DA in AT: SEXT Chro (fitted) ' num2str(c)]); %,'HandleVisibility','off'
plot(XX,ZZ,'b*-','MarkerSize',10,'LineWidth',2,'DisplayName', 'DA in AT: SEXT + QUAD Fringing Field');
%plot(XX,ZZ,'bo-','DisplayName', 'DA in AT SEXT'); %,'HandleVisibility','off'
%plot(DA_beta(:,1),DA_beta(:,2),'bo-','MarkerSize',7,'LineWidth',2,'DisplayName', 'DA in BETA: SEXT');
%plot(DA_beta_quadFF(:,1),DA_beta_quadFF(:,2),'m*-','MarkerSize',10,'LineWidth',2,'DisplayName', 'DA in BETA: SEXT + QUAD Fringing Field');
%plot(DA_beta_multip(:,1),DA_beta_multip(:,2),'gd-','MarkerSize',10,'LineWidth',2,'DisplayName', 'DA in BETA: MULTIPOLES');
plot(DA_beta_multip_Qff(:,1),DA_beta_multip_Qff(:,2),'rd-','MarkerSize',10,'LineWidth',2,'DisplayName', 'DA in BETA: MULTIPOLES + QUAD Fringing Field');
plot(XX_multip,ZZ_multip,'c*-','MarkerSize',9,'LineWidth',2,'DisplayName', ['DA in AT: QUAD/SEXT MULTIPOELS + QUAD Fringing Field ' ]); %,'HandleVisibility','off'
%plot(XX_quadFF,ZZ_quadFF,'b*-','DisplayName', 'DA in AT SEXT QUADff');
%plot(XX0,ZZ0,'b*-','DisplayName', 'DA in AT Chro 0');
% plot(XXq,ZZq,'mo-');
% plot(XXq0,ZZq0,'m*-');
xlabel('x [m]')
ylabel('z [m]')
set(gca,'FontSize',20)
set(gcf,'color','w')
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
xlim([-0.04 0.04])
ylim([0 0.025])
print('DA_Multipoles','-dpng','-r300')

%%

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
plot(x,y,'k-','DisplayName', 'Scaled vacuum chamber')
hold on; 
plot(x_bpipe,y_bpipe,'k--','DisplayName', 'Vacuum chamber')
plot(XX,ZZ,'g*-','MarkerSize',9,'LineWidth',2,'DisplayName', 'DA 3.17/1.64 Chro 0/0'); %,'HandleVisibility','off'
%plot(XX2,ZZ2,'c*-','MarkerSize',9,'LineWidth',2,'DisplayName', 'DA 3.16/1.58 Chro 2/2'); %,'HandleVisibility','off'
plot(XX_quadFF,ZZ_quadFF,'mo-','MarkerSize',9,'LineWidth',2,'DisplayName', 'DA 3.17/1.64 Chro 0/0: QUAD Fringing Field'); %,'HandleVisibility','off'
%plot(XX2_quadFF,ZZ2_quadFF,'yo-','MarkerSize',9,'LineWidth',2,'DisplayName', 'DA 3.16/1.58 Chro 2/2: QUAD Fringing Field'); %,'HandleVisibility','off'
%plot(DA_beta(:,1),DA_beta(:,2),'bo-','MarkerSize',7,'LineWidth',2,'DisplayName', 'DA in BETA: SEXT');
%plot(DA_beta_quadFF(:,1),DA_beta_quadFF(:,2),'m*-','MarkerSize',10,'LineWidth',2,'DisplayName', 'DA in BETA: SEXT + QUAD Fringing Field');
%plot(DA_beta_multip(:,1),DA_beta_multip(:,2),'gd-','MarkerSize',10,'LineWidth',2,'DisplayName', 'DA in BETA: MULTIPOLES');
%plot(DA_beta_multip_Qff(:,1),DA_beta_multip_Qff(:,2),'rd-','MarkerSize',10,'LineWidth',2,'DisplayName', 'DA in BETA: MULTIPOLES + QUAD Fringing Field');
plot(XX_multip_quadFF,ZZ_multip_quadFF,'bd-','MarkerSize',9,'LineWidth',2,'DisplayName', ['DA 3.17/1.64 Chro 0/0: QUAD/SEXT MULTIPOELS + QUAD Fringing Field ' ]); %,'HandleVisibility','off'
%plot(XX2_multip_quadFF,ZZ2_multip_quadFF,'rd-','MarkerSize',9,'LineWidth',2,'DisplayName', ['DA 3.16/1.58 Chro 2/2: QUAD/SEXT MULTIPOELS + QUAD Fringing Field ' ]); %,'HandleVisibility','off'
%plot(XX_quadFF,ZZ_quadFF,'b*-','DisplayName', 'DA in AT SEXT QUADff');
%plot(XX0,ZZ0,'b*-','DisplayName', 'DA in AT Chro 0');
% plot(XXq,ZZq,'mo-');
% plot(XXq0,ZZq0,'m*-');
xlabel('x [m]')
ylabel('z [m]')
set(gca,'FontSize',20)
set(gcf,'color','w')
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
xlim([-0.04 0.04])
ylim([0 0.03])
%print('DA_Multipoles_2lattice','-dpng','-r300')

%%



