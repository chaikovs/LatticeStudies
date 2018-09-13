clear all
close all
%addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors')

% load lattice
%load ../../ESRFLattice.mat
% load lattice
ring = thomx_at2();

% get indexes
indm=find(atgetcells(ring,'Class','Monitor'));

% get indexes
indq=find(atgetcells(ring,'Class','Quadrupole'));
indd=find(atgetcells(ring,'Class','bend'));
% 
rerr=atsetrandomerrors(...
    ring,...
    indd,...
    indm,...
    1,...
    1e-4,...
    2.5,...
    'dpb1');

% rerr=atsetrandomerrors(...
%     ring,...
%     indq,...
%     indm,...
%     1,...
%     1e-2,...
%     2.5,...
%     'dpb2');

figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(rerr,'comment',[],@plPolynomBComp);
%saveas(gca,'FieldIntegral.fig')
%export_fig('FieldIntegral.jpg','-r300')


% plots
figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(ring)
%saveas(gca,'OrbitNoErr.fig')
%export_fig('OrbitNoErr.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(rerr)
%saveas(gca,'OrbitWithErr.fig')
%export_fig('OrbitWithErr.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(rerr,'comment',[],@plClosedOrbit)
%saveas(gca,'OrbitWithErr.fig')
%export_fig('OrbitWithErr.jpg','-r300')
