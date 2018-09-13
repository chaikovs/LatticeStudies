
%%
ring = ThomX_017_064_r56_02_chro00_AT2;
%%
indBPM=find(atgetcells(ring,'FamName','BPMx'));
indHCorT=find(atgetcells(ring,'Class','Corrector'));
indHCor = indHCorT(1:2:end);
indVCorT=find(atgetcells(ring,'Class','Corrector'));
indVCor = indVCorT(1:2:end);

sBPM=findspos(ring,indBPM);
sRING=findspos(ring,1:length(ring)+1);
%%

ringerr = atsetfieldvalues(ring, findcells(ring,'FamName','HCOR'), 'KickAngle',{1,1},0.1e-3);
ringerr1 = atsetfieldvalues(ring,68 , 'KickAngle',{1,1},0.1e-3);
%%

orbit4err = findorbit4(ringerr1, 0,indBPM); % 1:length(rcor)+1
xorbiterr=orbit4err(1,:);
yorbiterr=orbit4err(3,:);

figure(13);
plot(sBPM,1e3*xorbiterr,'b.-','MarkerSize',8)
hold on; 
plot(sBPM,1e3*yorbiterr,'r.-');
set(gcf,'color','w')
set(gca,'fontsize',16');
 xlim([0 18])
% ylim([0 6])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('s-position [m]');
ylabel('y [mm]');

%%

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(ringerr1,'comment',[],@plClosedOrbit)

