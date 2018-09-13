
clear all; close all
%%

ring = ThomX_FF();

%%

ring = thomx_at2();

%%

% ring = test();
% 
% ring=atsetfieldvalues(ring,findcells(ring,'PassMethod','CorrectorPass'),...
%     'PassMethod','StrMPoleSymplectic4Pass');
%%

indq=find(atgetcells(ring,'Class','Quadrupole'));

%indBPM=find(atgetcells(ring,'Class','Monitor'));
indBPM=find(atgetcells(ring,'FamName','BPMx'));

indHCor=find(atgetcells(ring,'FamName','HCOR'));
indVCor=find(atgetcells(ring,'FamName','VCOR'));


indQCor=find(atgetcells(ring,'FamName','HCOR'));

%indHCor=find(atgetcells(ring,'Class','Corrector'));
%indHCor=find(atgetcells(ring,'Class','Sextupole'));
%indVCor=find(atgetcells(ring,'Class','Corrector'));
%indVCor=find(atgetcells(ring,'Class','Sextupole'));

% indHCor=find(atgetcells(ring,'FamName','CHV0'));
% indVCor=find(atgetcells(ring,'FamName','CHV1'));
% indSCor=find(atgetcells(ring,'Class','Corrector'));
% indQCor=find(atgetcells(ring,'Class','Corrector'));

%%

%findcells(ring,'FamName','BEND')

%%

dx=10e-6*randn(size(indq));
dy=10e-6*randn(size(indq));

ringerr=atsetshift(ring,indq,dx,dy);

%%
% define s-axis rotation errors
dt=1e-6*randn(size(indq)); % random errors of 1um

% set errors
ringerr=atsettilt(ring,indq,dt);

%%

figure('units','normalized','position',[0.1 0.4 0.65 0.35])
atplot(ringerr,[0,10],'comment',[],@plPolynomBComp);


%%
figure
atplot(ring,'comment',[],@plClosedOrbit)
figure
atplot(ringerr,'comment',[],@pltmisalignments)
figure
atplot(ringerr,'comment',[],@plClosedOrbit)
%%

%compute response matrix if it doesn't already exist

if ~exist('ModelRM')
ModelRM...
        =getresponsematrices(...
        ring,...
        indBPM,...
        indHCor,...
        indVCor,...
        [],...
        indQCor,...
        []',...
        [0 0 0 0 0 0]',...
        [1 2 3]);
end

%%

%now orbit correction
[rcor,inCOD,hs,vs]=atcorrectorbit(ringerr,...
    indBPM,...
    indHCor',...
    indVCor',...
    [0 0 0 0 0 0]',...
    [10 10],...
    [false true],...
    1.0,...
    ModelRM,...
    zeros(2,length(indBPM)),...
    [],...
    true);

%%

%find new closed orbit
o=findorbit6Err(ringerr,indBPM,inCOD);
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
%saveas(gca,'OrbitCor.fig');
%export_fig('OrbitCor.jpg','-r300');


% plot output

figure;
subplot(2,1,1);bar(hs);ylabel('hor.')
subplot(2,1,2);bar(vs);ylabel('ver.')

inCOD(5)

rcor0=rcor;
%%

%r_new = atsetfieldvalues(ring,indHCor(1),'KickAngle',1e-4); %atsetfieldvalues
r_new = atsetfieldvalues(ring,indVCor(1),'KickAngle',{1,2},1e-5); %atsetfieldvalues

%%

atplot(r_new,'comment',[],@plClosedOrbit)

%%
