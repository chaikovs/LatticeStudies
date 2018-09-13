clear all
close all

% addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors')
% addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/random')
% addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/errordisplayfunctions');

% load lattice
ring = thomx_at2_BPM();

r0=ring;

% define errors to set
ie=1;
% sextupoles
inds=findcells(r0,'Class','Sextupole');
errstruct(ie).indx=inds;
errstruct(ie).type='psi'; % roll
errstruct(ie).sigma=200*1e-6;
ie=ie+1;

indqm=[findcells(r0,'Class','Quadrupole')];
errstruct(ie).indx=indqm;
errstruct(ie).type='x';
errstruct(ie).sigma=150*1e-6;
ie=ie+1;
errstruct(ie).indx=indqm;
errstruct(ie).type='y';
errstruct(ie).sigma=170*1e-6;
ie=ie+1;

%% set errors

magindex=arrayfun(@(a)a.indx,errstruct,'un',0);
type=arrayfun(@(a)a.type,errstruct,'un',0);
sigma=arrayfun(@(a)a.sigma,errstruct,'un',0);

rerr=atsetrandomerrors(...
    r0,...
    magindex,...
    findcells(r0,'Class','Monitor'),...
    123456,...
    sigma,...
    2.5,...
    type);


figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(rerr,[0,9],'comment',[],@pltmisalignments);
%saveas(gca,'LargeList.fig')
%export_fig('LargeList.jpg','-r300')

