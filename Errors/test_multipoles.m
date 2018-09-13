
clear all
close all

% load lattice
%load ../../ESRFLattice.mat
 %ring = ThomX_017_064_r56_02_chro00;
%ring = ThomX_017_058_r56_02_chro22;

ring = ThomX_017_064_r56_02_chro00_AT2;


% get indexes
indq=find(atgetcells(ring,'Class','Quadrupole'));
inds=find(atgetcells(ring,'Class','Sextupole'));

%QUAD
% b6   -20e-4
% b10  -23e-4   
% b14  -24e-4 

%SEXT
% b9    -3.29e-3
% b15   -1.16e-3

% set multipole errors
bn_quad=[0 0 0 0 0 -19.54 0 0 0 -23.46 0 0 0 -24.15]*1e-4; % 

bn_sext=[0 0 0 0 0 0 0 0 -3.29 0 0 0 0 0 -1.16]*1e-3; % 

% an=[
%     0
%     0
%        ]'*1e-4;

[rerr_quad,PolB_quad,PolA_quad]=AssignFieldErr(ring,indq,2,20*1e-3,bn_quad);
[rerr,PolB,PolA]=AssignFieldErr(rerr_quad,inds,3,20*1e-3,bn_sext);

%%

% atwritem(rerr, '/Users/ichaikov/Documents/MATLAB/thomx-mml/machine/THOMX/StorageRing/Lattices/ThomX_016_058_r56_02_chro22_multip')
%atwritem(rerr, '/Users/ichaikov/Documents/MATLAB/thomx-mml/machine/THOMX/StorageRing/Lattices/ThomX_017_064_r56_02_chro00_multip_AT2')

%%

DA_beta = load('/Users/ichaikov/Documents/MATLAB/thomx-mml/machine/THOMX/StorageRing/DA/DA_SXfit1');
[XX,ZZ]   = atdynap(ring, 500,0,0.02); 
[XX_mp,ZZ_mp]   = atdynap(rerr, 500,0,0.02); 

[l,t,c] = atlinopt(ring,0,1);
[l_mp,t_mp,c_mp] = atlinopt(rerr,0,1);

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


figure('units','normalized','position',[0.3 0.3 0.45 0.35])
plot(x,y,'k-','DisplayName', 'Scaled vacuum chamber')
hold on; 
plot(x_bpipe,y_bpipe,'k--','DisplayName', 'Vacuum chamber')
%plot(XX,ZZ,'bo-','MarkerSize',7,'LineWidth',2,'DisplayName', ['DA in AT SEXT Chro '  num2str(c)]); %,'HandleVisibility','off'
plot(DA_beta(:,1),DA_beta(:,2),'r*-','MarkerSize',10,'LineWidth',2,'DisplayName', 'DA in BETA: SEXT Chro 00');
plot(XX,ZZ,'bo-','MarkerSize',7,'DisplayName', 'DA in AT SEXT Chro 00'); %,'HandleVisibility','off'
plot(XX_mp,ZZ_mp,'mo-','MarkerSize',7,'LineWidth',2,'DisplayName', ['DA in AT: QUAD/SEXT MULTIPOLES Chro (fitted) ' num2str(c_mp)]); %,'HandleVisibility','off'
xlabel('x [m]')
ylabel('z [m]')
set(gca,'FontSize',20)
set(gcf,'color','w')
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
grid on
xlim([-0.05 0.05])
ylim([0 0.025])
addlabel(1, 0, datestr(clock,0))
%print('DA_BETA_AT_chro00','-dpng','-r300')


