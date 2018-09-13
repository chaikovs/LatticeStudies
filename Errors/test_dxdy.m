clear all
close all


% load lattice
%load ../../ESRFLattice.mat
% load lattice
ring = ThomX_017_064_r56_02_chro00_AT2();

% get indexes
indm=find(atgetcells(ring,'FamName','BPMx'));

% get indexes
indq=find(atgetcells(ring,'Class','Quadrupole'));
indd=find(atgetcells(ring,'Class','bend'));
%indd=find(atgetcells(ring,'FamName','SX2'));
%indsext=[findcells(r0,'Class','Sextupoles')];

% define alignemnt errors
% dx=100e-6*randn(size(indq)); % random errors of 10um
% dy=100e-6*randn(size(indq));

dx=100e-6*randn(size(indd)); % random errors of 10um
dy=100e-6*randn(size(indd));

% set errors
% ringerr=atsetshift(ring,indq,dx,dy);
ringerr=atsetshift(ring,indd,dx,dy);

% plots
figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(ring,'comment',[],@plClosedOrbit)
%saveas(gca,'OrbitNoErr.fig')
%export_fig('OrbitNoErr.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(ringerr,[0, 9],'comment',[],@plClosedOrbit)
%saveas(gca,'OrbitWithErr.fig')
%export_fig('OrbitWithErr.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(ringerr,'comment',[],@pltmisalignments)
%saveas(gca,'SetErrDxDy.fig')
%export_fig('SetErrDxDy.jpg','-r300')

return
%%
figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(ringerr,'comment',[],@plot_betabeat)