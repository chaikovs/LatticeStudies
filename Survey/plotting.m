global GLOBVAL 
GLOBVAL.E0=50E6;
GLOBVAL.LatticeFile='test';

RING = TDR_017_064_r56_02_Dff412_chro00;
 

%%

figure('units','normalized','position',[0.1 0.4 0.65 0.25])
set(gcf,'color','w')
drawlattice
%set(h2,'YTick',[])
set(gca,'FontSize',18)
xlim([0 9])
%xlabel('s - position [m]');
%linkaxes([h1 h2],'x')
%set([h1 h2],'XGrid','On','YGrid','On');
ax = gca
ax.Visible = 'off'
% The axes is removed from PDF file, too.
%print(gcf,'-dpdf','-r300');
%print('thomxSR_lat.png','-dpng','-r300')