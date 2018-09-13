%error and correction for HEPS
ring = heps_Cor_BPM(1);

% get indices
indq=find(atgetcells(ring,'Class','Quadrupole'));
indBPM=find(atgetcells(ring,'FamName','BPMx'));
%indBPM=find(atgetcells(ring,'Class','Monitor'));

indHCor=find(atgetcells(ring,'Class','Corrector'));
%indHCor=find(atgetcells(ring,'Class','Sextupole'));
indVCor=find(atgetcells(ring,'Class','Corrector'));
%indVCor=find(atgetcells(ring,'Class','Sextupole'));


% indHCor=find(atgetcells(ring,'FamName','CHV0'));
% indVCor=find(atgetcells(ring,'FamName','CHV1'));
indSCor=find(atgetcells(ring,'Class','Corrector'));
indQCor=find(atgetcells(ring,'Class','Corrector'));

% define alignemnt errors
dx=1e-6*randn(size(indq)); % random errors of 1um
dy=1e-6*randn(size(indq));

% set errors
rerr=atsetshift(ring,indq,dx,dy);


%compute response matrix if it doesn't already exist

if ~exist('ModelRM')
ModelRM...
        =getresponsematrices(...
        ring,...
        indBPM,...
        indHCor,...
        indVCor,...
        indSCor,...
        indQCor,...
        [],...
        [0 0 0 0 0 0]',...
        [1 2 3]);
end
    
%now orbit correction
[rcor,inCOD,hs,vs]=atcorrectorbit(rerr,...
    indBPM,...
    indHCor',...
    indVCor',...
    [0 0 0 0 0 0]',...
    [50 50],...
    [false true],...
    1.0,...
    ModelRM,...
    zeros(2,length(indBPM)),...
    [],...
    true);


%find new closed orbit
o=findorbit6Err(rerr,indBPM,inCOD);
oxe=o(1,:);
oye=o(3,:);

o=findorbit6Err(rcor,indBPM,inCOD);
oxc=o(1,:);
oyc=o(3,:);

sBPM=findspos(rcor,indBPM);
figure;subplot(2,1,1);
plot(sBPM,oxe,'.-');hold on; plot(sBPM,oxc,'.-');
legend('before','after');
xlabel('s [m]');
ylabel('hor. COD');
subplot(2,1,2);
plot(sBPM,oye,'.-');hold on; plot(sBPM,oyc,'.-');
legend('before','after');
xlabel('s [m]');
ylabel('ver. COD');
saveas(gca,'OrbitCor.fig');
export_fig('OrbitCor.jpg','-r300');


% plot output

figure;
subplot(2,1,1);bar(hs);ylabel('hor.')
subplot(2,1,2);bar(vs);ylabel('ver.')

inCOD(5)

rcor0=rcor;

%%
