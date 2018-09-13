clear all
close all

% load lattice
ring = thomx_at2_BPM();

% get indexes
indHCor=find(atgetcells(ring,'FamName','HCOR'));
indc=indHCor;

%set correctors to 1e-5 rad horizontal kick
ring0=ring;

ring=atsetfieldvalues(ring,indc(3),'PolynomB',{1,1},1e-2);


% define s-axis rotation errors in correctors
dt=1e-6*randn(size(indc)); % random errors of 1um

% set errors
ringerr=atsettilt(ring,indc,dt);


% plots
[twin,~,~]=atlinopt(ring0,0,1);

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(ring,'comment',[],'inputtwiss',twin,@plClosedOrbit)
%saveas(gca,'OrbitNoErrTiltCor.fig')
%export_fig('OrbitNoErrTiltCor.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(ringerr,'comment',[],'inputtwiss',twin,@plClosedOrbit)
%saveas(gca,'OrbitWithErrTiltCor.fig')
%export_fig('OrbitWithErrTiltCor.jpg','-r300')

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
s=findspos(ring,indc);
plot(s,dt,'b.');
atplotsyn(gca,ringerr);
axis tight;
xlabel('s [m]');
ylabel('s-axis rotation [rad]');

%saveas(gca,'SetErrTiltCor.fig')
%export_fig('SetErrTiltCor.jpg','-r300')

