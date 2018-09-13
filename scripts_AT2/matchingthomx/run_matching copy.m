close all; clear all

r_original=thomx04_original();
r_ml = thomxML();

%%

r_optic=matchthomx_diplength(r_ml);


%%

figure; c=atplot(r_optic);   % after match
%%

k_original = getkval(r_original);
k_match = getkval(r_optic);
%%

% plot 
% figure; c=atplot(r_original);  c.comment.String{3}=[' alpha= ' num2str(mcf(r))];% before match
% figure; c=atplot(r);  c.comment.String{3}=[' alpha= ' num2str(mcf(r_optic))]; % after match

figure; c=atplot(r_ml);  %print('beta_atplot_ML.png', '-dpng', '-r300');% before match
figure; c=atplot(r_optic); %print('beta_atplot_Match.png', '-dpng', '-r300');  % after match
figure; c=atplot(r_original); %print('beta_atplot_orig.png', '-dpng', '-r300');  % after match
%%

figure
subplot 211
set(gca,'FontSize',18)
plot(k_match(1,:),'ro-','DisplayName', 'After matching')
hold on
plot(k_original(1,:),'bo-','DisplayName', 'Before matching')
hold off
legend('show','Location','NorthWest');
xlabel('Family number');
ylabel('K value [1/m^2]');
subplot 212
set(gca,'FontSize',18)
bar(k_match(1,:)./k_original(1,:))
xlabel('Family number');
ylabel('K_{match}/K_{orig}');
print('K_matching.png', '-dpng', '-r300');
%%
%indBPM=find(atgetcells(r_optic,'FamName','BPMx'))';
indBPM=1:length(r_optic);

[l,t,c]=atlinopt(r_optic,0,indBPM);
[le,te,ce]=atlinopt(r_original,0,indBPM);
[lc,tc,cc]=atlinopt(r_ml,0,indBPM);
bx0=arrayfun(@(a)a.beta(1),l);
bxe=arrayfun(@(a)a.beta(1),le);
bxc=arrayfun(@(a)a.beta(1),lc);
figure;
plot(bx0); hold on; 
plot(bxe); 
plot(bxc);
legend('match','original','magn_length')

%%
% TDR_good_017_064_r56_02_DIP_MagnLength
% gettune
% 
% [BetaX, BetaY, BPMs] = modeltwiss('Beta');
% 
% figure
% h1 = subplot(5,1,[1 2]);
% plot(BPMs,BetaX,'.-b', 'Markersize',10);
% hold on
% xlim([0 BPMs(end)]);
% ylabel('\beta_x [m]');
% %title('\beta-functions');
% 
% h2 = subplot(5,1,3);
% drawlattice 
% set(h2,'YTick',[])
% 
% h3 = subplot(5,1,[4 5]);
% plot(BPMs,BetaY,'.-r', 'Markersize',10);
% hold on
% xlabel('s - position [m]');
% ylabel('\beta_z [m]');
% 
% linkaxes([h1 h2 h3],'x')
% set([h1 h2 h3],'XGrid','On','YGrid','On');
% 
% TDR_good_017_064_r56_02_DIP_MagnLength_Match
% gettune
% 
% [BetaX, BetaY, BPMs] = modeltwiss('Beta');
% 
% h1  = subplot(5,1,[1 2]);
% set(gca,'FontSize',16)
% plot(BPMs, BetaX, 'k.-', 'Markersize',10);
% hold off
% xlabel('s (m)'); ylabel('\beta_x [meters]')
% u = legend('Before matching','After matching');
% set(u,'Location','NorthEast','FontSize',12)
% 
% h2 = subplot(5,1,[4 5]);
% set(gca,'FontSize',16)
% plot(BPMs, BetaY, 'k.-', 'Markersize',10);
% xlabel('s (m)'); ylabel('\beta_z [meters]')
% u = legend('Before matching','After matching');
% set(u,'Location','NorthWest','FontSize',12)
% 
% %linkaxes([a1,a2],'x')
% hold off

%% Beta function


TDR_good_017_064_r56_04_sx_Dff_BPM_ML
gettune

[BetaX, BetaY, BPMs] = modeltwiss('Beta');

figure
h1 = subplot(5,1,[1 2]);
plot(BPMs,BetaX,'.-g', 'Markersize',10);
hold on
xlim([0 BPMs(end)]);
ylabel('\beta_x [m]');
%title('\beta-functions');

h2 = subplot(5,1,3);
drawlattice 
set(h2,'YTick',[])

h3 = subplot(5,1,[4 5]);
plot(BPMs,BetaY,'.-g', 'Markersize',10);
hold on
xlabel('s - position [m]');
ylabel('\beta_z [m]');

linkaxes([h1 h2 h3],'x')
set([h1 h2 h3],'XGrid','On','YGrid','On');

TDR_good_017_064_r56_02_DIP_MagnLength_Match
gettune

[BetaX, BetaY, BPMs] = modeltwiss('Beta');

h1 = subplot(5,1,[1 2]);
plot(BPMs,BetaX,'.-k', 'Markersize',10);

xlim([0 BPMs(end)]);
ylabel('\beta_x [m]');
%title('\beta-functions');

h2 = subplot(5,1,3);
drawlattice 
set(h2,'YTick',[])

h3 = subplot(5,1,[4 5]);
plot(BPMs,BetaY,'.-k', 'Markersize',10);

xlabel('s - position [m]');
ylabel('\beta_z [m]');

linkaxes([h1 h2 h3],'x')
set([h1 h2 h3],'XGrid','On','YGrid','On');


TDR_good_017_064_r56_04_sx_Dff_BPM
gettune

[BetaX, BetaY, BPMs] = modeltwiss('Beta');

h1  = subplot(5,1,[1 2]);
set(gca,'FontSize',16)
plot(BPMs, BetaX, 'b.-', 'Markersize',10);
hold off
xlabel('s (m)'); ylabel('\beta_x [meters]')
u = legend('Before matching','After matching','Original');
set(u,'Location','NorthEast','Orientation','horizontal','FontSize',12)

h2 = subplot(5,1,[4 5]);
set(gca,'FontSize',16)
plot(BPMs, BetaY, 'r.-', 'Markersize',10);
xlabel('s (m)'); ylabel('\beta_z [meters]')
u = legend('Before matching','After matching','Original');
set(u,'Location','NorthWest','Orientation','horizontal','FontSize',12)

%linkaxes([a1,a2],'x')
hold off
print('beta_matching.png', '-dpng', '-r300');

%% Dispersion

TDR_good_017_064_r56_04_sx_Dff_BPM
[Dx, Dy, Sx, Sy] = modeldisp

% Figure to check
% plot betax and betay in two subplots
figure(1)
h1 = subplot(5,1,[1 4]);
set(gca,'FontSize',16)
plot(Sx,Dx,'.-b', 'Markersize',10, 'Linewidth', 1.6)
hold on
xlim([0 Sx(end)]);
ylabel('\eta_x [m]');
%title('Optical-functions');
h2 = subplot(5,1,5);
set(gca,'FontSize',16)
drawlattice 
set(h2,'YTick',[])
xlabel('s - position [m]');

linkaxes([h1 h2],'x')
set([h1 h2],'XGrid','On','YGrid','On');


TDR_good_017_064_r56_04_sx_Dff_BPM_ML

h1 = subplot(5,1,[1 4]);
set(gca,'FontSize',16)
plot(Sx, Dx,'.-g', 'Markersize',10, 'Linewidth', 1.6)
xlim([0 Sx(end)]);
ylabel('\eta_x [m]');
%title('Optical-functions');
h2 = subplot(5,1,5);
set(gca,'FontSize',16)
drawlattice 
set(h2,'YTick',[])
xlabel('s - position [m]');

linkaxes([h1 h2],'x')
set([h1 h2],'XGrid','On','YGrid','On');

TDR_good_017_064_r56_04_sx_Dff_BPM_Match

h1 = subplot(5,1,[1 4]);
set(gca,'FontSize',16)
plot(Sx, Dx,'.-k', 'Markersize',10, 'Linewidth', 1.6)
xlim([0 Sx(end)]);
ylabel('\eta_x [m]');
%title('Optical-functions');
u = legend({'Original \eta_x','Magn Length \eta_x','Match \eta_x'});
set(u,'Location','NorthEast','Orientation','horizontal','FontSize',12)
h2 = subplot(5,1,5);
set(gca,'FontSize',16)
drawlattice 
set(h2,'YTick',[])
xlabel('s - position [m]');

linkaxes([h1 h2],'x')
set([h1 h2],'XGrid','On','YGrid','On');
