clear all
close all

addpath('/Users/ichaikov/Documents/MATLAB/thomx-mml/at_v2.0/atmat/pubtools/LatticeTuningFunctions/errors');
addpath('/Users/ichaikov/Documents/MATLAB/thomx-mml/at_v2.0/atmat/pubtools/LatticeTuningFunctions/errors/random')

% load lattice
ring = ThomX_017_064_r56_02_chro00_AT2();

% get indexes
indm=find(atgetcells(ring,'FamName','BPMx'));

indq=find(atgetcells(ring,'Class','Quadrupole')); % girders are defined by GS and GE markers (start and end of girder)

rerr=atsetrandomerrors(...
    ring,...
    indq,...
    indm,...
    1,...
    1e-4,...
    2.5,...
    's');

figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(ring)
hold on;
atplotsyn(gca,rerr);

%saveas(gca,'DeltaSQuad.fig')
%export_fig('DeltaSQuad.jpg','-r300')
figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(rerr)

%% get indexes
%indm=find(atgetcells(ring,'Class','Monitor'));
indm=find(atgetcells(ring,'FamName','BPMx'));
indd=find(atgetcells(ring,'FamName','BEND'));
%indd=find(atgetcells(ring,'Class','Bend')); % girders are defined by GS and GE markers (start and end of girder)

rerr=atsetrandomerrors(...
    ring,...
    indd(1:1:end),...
    indm,...
    123456,...
    1e-2,...
    2.5,...
    's');

% rerr=atsetrandomerrors(...
%     ring,...
%     indd(1:1:end),...
%     indm,...
%     123456,...
%     1e-4,...
%     2.5,...
%     'dpb1');

figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(ring)
hold on;
atplotsyn(gca,rerr);

%saveas(gca,'DeltaSDipZoom.fig')
%export_fig('DeltaSDipZoom.jpg','-r300')

figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(rerr)


figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(rerr,'comment',[],@plClosedOrbit)
