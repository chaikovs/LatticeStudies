
%%
% global THERING
% indx=1:length(THERING);  
% T=twissring(THERING,0,indx);   
% beta=cat(1,T.beta);
%[bx by] = modelbeta;
%%

%ring0= ThomX_017_064_r56_02(); 
ring = ThomX_017_064_r56_02_chro00();
ring_multip = ThomX_017_064_r56_02_chro00_multip_AT2();
%ring_multip = ThomX_017_064_r56_02_chro00_multip();
% ring_chro00= ThomX_017_064_r56_02_chro00(); 
%ring2= ThomX_016_058_r56_02_chro22(); 
%ring2_multip= ThomX_016_058_r56_02_chro22_multip(); 

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

dp = 0.01;
[XX,ZZ]   = atdynap(ring, 1000,0,0.02); 
[XXplus,ZZplus]   = atdynap(ring, 1000,dp,0.02); 
[XXminus,ZZminus]   = atdynap(ring, 1000,-dp,0.02); 
%[XX_multip,ZZ_multip]   = atdynap(ring_multip, 500,0,0.02); 

% [XX2,ZZ2]   = atdynap(ring2, 500,0,0.02); 
% [XX2_multip,ZZ2_multip]   = atdynap(ring2_multip, 500,0,0.02);

[l,t,c] = atlinopt(ring,0,1);
%[l_mp,t_mp,c_mp] = atlinopt(ring_multip,0,1);

%%

% ring_quadFF =atsetfieldvalues(ring,find(atgetcells(ring,'Class','Quadrupole')),...
%     'PassMethod','QuadMPoleFringePass' );

%ring_multip_quadFF =atsetfieldvalues(ring_multip,find(atgetcells(ring_multip,'Class','Quadrupole')),...
%    'PassMethod','QuadMPoleFringePass' );

%ring2_quadFF =atsetfieldvalues(ring2,find(atgetcells(ring2,'Class','Quadrupole')),...
%    'PassMethod','QuadMPoleFringePass' );

%ring2_multip_quadFF =atsetfieldvalues(ring2_multip,find(atgetcells(ring2_multip,'Class','Quadrupole')),...
%    'PassMethod','QuadMPoleFringePass' );

% [XX_quadFF,ZZ_quadFF]   = atdynap(ring_quadFF, 500,0.0,0.02); 
% [XXplus_quadFF,ZZplus_quadFF]   = atdynap(ring_quadFF, 500,dp,0.02);
% [XXminus_quadFF,ZZminus_quadFF]   = atdynap(ring_quadFF, 500,-dp,0.02);
%[XX_multip_quadFF,ZZ_multip_quadFF]   = atdynap(ring_multip_quadFF, 500,0.0,0.02); 

%[XX2_quadFF,ZZ2_quadFF]   = atdynap(ring2_quadFF, 500,0.0,0.02); 
%[XX2_multip_quadFF,ZZ2_multip_quadFF]   = atdynap(ring2_multip_quadFF, 500,0.0,0.02); 

%DA_beta_quadFF = load('DA_SXfit1_Qff');

% [l_quadFF,t_quadFF,c_quadFF] = atlinopt(ring_quadFF,0,1);


[XX_multip,ZZ_multip]   = atdynap(ring_multip, 500,0.0,0.02); 
[XXplus_multip,ZZplus_multip]   = atdynap(ring_multip, 500,dp,0.02);
[XXminus_multip,ZZminus_multip]   = atdynap(ring_multip, 500,-dp,0.02);

[l_multip,t_multip,c_multip] = atlinopt(ring_multip,0,1);

%%
sigma_x = sqrt(8e-6/50/0.511*3.5);
sigma_z = sqrt(8e-6/50/0.511*2.7);

figure('units','normalized','position',[0.3 0.3 0.4 0.35])
plot(x,y,'k-','DisplayName', 'Scaled vacuum chamber')
hold on; 
plot(x_bpipe,y_bpipe,'k--','DisplayName', 'Vacuum chamber')
plot(XX,ZZ,'b.-','MarkerSize',12,'LineWidth',3,'DisplayName', ['DA in AT SEXT Chro '  num2str(c)]); %,'HandleVisibility','off'
plot(XXminus,ZZminus,'r.--','MarkerSize',12,'LineWidth',3,'DisplayName', 'SEXT, dp/p = - 1%');
%plot(XX0,ZZ0,'mo-','MarkerSize',7,'DisplayName', 'DA in AT SEXT Chro 00'); %,'HandleVisibility','off'
plot(XXplus,ZZplus,'m.--','MarkerSize',12,'LineWidth',3,'DisplayName', 'SEXT, dp/p = +1%'); %,'HandleVisibility','off'
xlabel('x [m]')
ylabel('z [m]')
set(gca,'FontSize',20)
set(gcf,'color','w')
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
xlim([-0.03 0.03])
ylim([0 0.02])
addlabel(1, 0, datestr(clock,0))
print('DA_SEXT_AT_chro00_dpp','-dpng','-r300')

figure('units','normalized','position',[0.3 0.3 0.4 0.35])
plot(x./sigma_x,y./sigma_z,'k-','DisplayName', 'Scaled vacuum chamber')
hold on; 
plot(x_bpipe./sigma_x,y_bpipe./sigma_z,'k--','DisplayName', 'Vacuum chamber')
plot(XX./sigma_x,ZZ./sigma_z,'b.-','MarkerSize',12,'LineWidth',3,'DisplayName', ['DA in AT SEXT Chro '  num2str(c)]); %,'HandleVisibility','off'
plot(XXminus./sigma_x,ZZminus./sigma_z,'r.--','MarkerSize',12,'LineWidth',3,'DisplayName', 'SEXT, dp/p = - 1%');
%plot(XX0,ZZ0,'mo-','MarkerSize',7,'DisplayName', 'DA in AT SEXT Chro 00'); %,'HandleVisibility','off'
plot(XXplus./sigma_x,ZZplus./sigma_z,'m.--','MarkerSize',12,'LineWidth',3,'DisplayName', 'SEXT, dp/p = +1%'); %,'HandleVisibility','off'
xlabel('x/\sigma_x')
ylabel('z/\sigma_z')
set(gca,'FontSize',20)
set(gcf,'color','w')
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
 xlim([-30 30])
 ylim([0 20])
addlabel(1, 0, datestr(clock,0))
print('DA_SEXT_AT_chro00_dpp_sigma','-dpng','-r300')
%%


figure('units','normalized','position',[0.3 0.3 0.4 0.35])
plot(x./sigma_x,y./sigma_z,'k-','DisplayName', 'Scaled vacuum chamber')
hold on; 
plot(x_bpipe./sigma_x,y_bpipe./sigma_z,'k--','DisplayName', 'Vacuum chamber')
%plot(XX,ZZ,'bo-','MarkerSize',7,'DisplayName', 'DA in AT SEXT'); %,'HandleVisibility','off'
%plot(DA_beta(:,1),DA_beta(:,2),'ro-','MarkerSize',7,'DisplayName', 'DA in BETA SEXT');
plot(XX./sigma_x,ZZ./sigma_z,'b.-','MarkerSize',10,'LineWidth',3,'DisplayName', ['SEXT + QUAD Fringing Field']); %,'HandleVisibility','off'
%plot(DA_beta_quadFF(:,1),DA_beta_quadFF(:,2),'r*-','MarkerSize',10,'LineWidth',2,'DisplayName', 'DA in BETA: SEXT + QUAD Fringing Field');
plot(XX_multip./sigma_x,ZZ_multip./sigma_z,'g-','MarkerSize',10,'LineWidth',3,'DisplayName', 'QUAD/SEXT MULTIPOELS + QUAD Fringing Field');
plot(XXplus_multip./sigma_x,ZZplus_multip./sigma_z,'m.--','MarkerSize',12,'LineWidth',3,'DisplayName', 'QUAD/SEXT MULTIPOELS + QUAD Fringing Field, dp/p = +1%');
plot(XXminus_multip./sigma_x,ZZminus_multip./sigma_z,'r.--','MarkerSize',12,'LineWidth',3,'DisplayName', 'QUAD/SEXT MULTIPOELS + QUAD Fringing Field, dp/p = -1%');
%plot(XX0,ZZ0,'b*-','DisplayName', 'DA in AT Chro 0');
% plot(XXq,ZZq,'mo-');
% plot(XXq0,ZZq0,'m*-');
xlabel('x/\sigma_x')
ylabel('z/\sigma_z')
set(gca,'FontSize',20)
set(gcf,'color','w')
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
xlim([-30 30])
 ylim([0 25])
addlabel(1, 0, datestr(clock,0))
print('DA_SEXT_AT_chro00_QuadFF_dpp_sigma','-dpng','-r300')

figure('units','normalized','position',[0.3 0.3 0.4 0.35])
plot(x,y,'k-','DisplayName', 'Scaled vacuum chamber')
hold on; 
plot(x_bpipe,y_bpipe,'k--','DisplayName', 'Vacuum chamber')
%plot(XX,ZZ,'bo-','MarkerSize',7,'DisplayName', 'DA in AT SEXT'); %,'HandleVisibility','off'
%plot(DA_beta(:,1),DA_beta(:,2),'ro-','MarkerSize',7,'DisplayName', 'DA in BETA SEXT');
plot(XX,ZZ,'b-','MarkerSize',12,'LineWidth',3,'DisplayName', ['SEXT + QUAD Fringing Field']); %,'HandleVisibility','off'
%plot(DA_beta_quadFF(:,1),DA_beta_quadFF(:,2),'r*-','MarkerSize',10,'LineWidth',2,'DisplayName', 'DA in BETA: SEXT + QUAD Fringing Field');
plot(XX_multip,ZZ_multip,'g-','MarkerSize',12,'LineWidth',3,'DisplayName', 'QUAD/SEXT MULTIPOELS + QUAD Fringing Field');
plot(XXplus_multip,ZZplus_multip,'m.--','MarkerSize',12,'LineWidth',3,'DisplayName', 'QUAD/SEXT MULTIPOELS + QUAD Fringing Field, dp/p = +1%');
plot(XXminus_multip,ZZminus_multip,'r.--','MarkerSize',12,'LineWidth',3,'DisplayName', 'QUAD/SEXT MULTIPOELS + QUAD Fringing Field, dp/p = -1%');
%plot(XX0,ZZ0,'b*-','DisplayName', 'DA in AT Chro 0');
% plot(XXq,ZZq,'mo-');
% plot(XXq0,ZZq0,'m*-');
xlabel('x [m]')
ylabel('z [m]')
set(gca,'FontSize',20)
set(gcf,'color','w')
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
xlim([-0.03 0.03])
ylim([0 0.02])
addlabel(1, 0, datestr(clock,0))
print('DA_SEXT_AT_chro00_QuadFF_dpp','-dpng','-r300')

%%






