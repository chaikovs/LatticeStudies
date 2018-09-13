
close all
clear all

%addpath(genpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/correction/'));
%addpath(genpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/'));

% load lattice
ring=THOMX();

indBPM=1:length(ring);   %find(atgetcells(ring,'FamName','BPMx'))';
s=findspos(ring,indBPM);

rerr=atsetcavity(ring,300e3,0,30); % beta =1 ;

o=findorbit6(ring,indBPM);

f0 = ring{2}.Frequency;

ring=atsetcavity_TESTBETANOT1(ring,300e3,0,30); 

o1=findorbit6(ring,indBPM);
f1 = ring{2}.Frequency;

figure;
set(gca,'FontSize',18)
plot(s,o(1,:));
hold on;
plot(s,o1(1,:));
legend(['beta=1, f_{RF}=' num2str(f0) ' Hz'],['beta thomx, f_{RF}=' num2str(f1) ' Hz'])
ylabel('x-orbit [m] findorbit6');
xlabel('s [m]');
title(['RF frequency difference ' num2str((f0-f1)/1e3) ' kHz'])
print('freqRF_diff.png', '-dpng', '-r300');

