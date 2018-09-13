

%% To test one error
% 
% ring = ThomX_017_064_r56_02_chro00_AT2();
% 
% %indm=find(atgetcells(ring,'Class','Monitor'));
% indm=find(atgetcells(ring,'FamName','BPMx'));
% indd=find(atgetcells(ring,'FamName','BEND'));
% % rerr=atsetrandomerrors(...
% %     ring,...
% %     indd(1:1:end),...
% %     indm,...
% %     123456,...
% %     1e-2,...
% %     2.5,...
% %     's');
% 
% rerr=atsetrandomerrors(...
%     ring,...
%     indd(1:1:end),...
%     indm,...
%     123456,...
%     100e-6,...
%     2.5,...
%     'x');
% 
% 
% figure('units','normalized','position',[0.3 0.3 0.45 0.35])
% atplot(rerr,'comment',[],@plClosedOrbit)

%%

clear all; close all;


% dx_dipole = 30e-6; %100*1e-6;
% dy_dipole = 30e-6;%100*1e-6;
% ds_dipole = 30e-6;%100*1e-6;
% tilt_dipole = 200e-6;%500*1e-6;
% dkk_dipole = 1e-3;%5e-3;%1.2e-3;
% 
% dx_quad = 30e-6;%100*1e-6;
% dy_quad = 30e-6;%100*1e-6;
% ds_quad = 30e-6;%100*1e-6;
% tilt_quad = 200e-6;%500*1e-6;
% dkk_quad= 1e-3;%5e-3;
% 
% dx_sext = 30e-6;%100*1e-6;
% dy_sext = 30e-6;%100*1e-6;
% ds_sext = 30e-6;%100*1e-6;
% tilt_sext = 200e-6;%500*1e-6;
% dkk_sext= 1e-3;%5e-3;

dx_dipole = 100*1e-6;
dy_dipole = 100*1e-6;
ds_dipole = 100*1e-6;
tilt_dipole = 500*1e-6;
dkk_dipole = 5e-3;%1.2e-3;

dx_quad =  100*1e-6;%100*1e-6;
dy_quad =  100*1e-6;%100*1e-6;
ds_quad =  100*1e-6;%100*1e-6;
tilt_quad = 500*1e-6;%500*1e-6;
dkk_quad= 5e-3;

dx_sext = 100*1e-6;
dy_sext = 100*1e-6;
ds_sext = 100*1e-6;
tilt_sext = 500*1e-6;
dkk_sext= 5e-3;


% load lattice
ring_woErr = ThomX_017_064_r56_02_chro00_AT2();
ring = ThomX_017_064_r56_02_chro00_multip_AT2();

r0=ring;

%QUAD FF
% ring_quadFF =atsetfieldvalues(r0,find(atgetcells(r0,'Class','Quadrupole')),...
%     'PassMethod','QuadMPoleFringePass' );
% r0 = ring_quadFF;

% define errors to set
ie=1;

%dipoles
inds=findcells(r0,'Class','bend');
errstruct(ie).indx=inds;
errstruct(ie).type='psi'; % roll
errstruct(ie).sigma=tilt_dipole;
ie=ie+1;
errstruct(ie).indx=inds;
errstruct(ie).type='x'; 
errstruct(ie).sigma=dx_dipole;
ie=ie+1;
errstruct(ie).indx=inds;
errstruct(ie).type='y'; 
errstruct(ie).sigma=dy_dipole;
ie=ie+1;
errstruct(ie).indx=inds;
errstruct(ie).type='s'; 
errstruct(ie).sigma=ds_dipole;
ie=ie+1;

errstruct(ie).indx=inds;
errstruct(ie).type='dpb1'; 
errstruct(ie).sigma=dkk_dipole;
ie=ie+1;

% Quadrupoles
indqm=[findcells(r0,'Class','Quadrupole')];
errstruct(ie).indx=indqm;
errstruct(ie).type='psi'; % roll
errstruct(ie).sigma=tilt_quad;
ie=ie+1;
errstruct(ie).indx=indqm;
errstruct(ie).type='x';
errstruct(ie).sigma=dx_quad;
ie=ie+1;
errstruct(ie).indx=indqm;
errstruct(ie).type='y';
errstruct(ie).sigma=dy_quad;
ie=ie+1;
errstruct(ie).indx=indqm;
errstruct(ie).type='s';
errstruct(ie).sigma=ds_quad;
ie=ie+1;

errstruct(ie).indx=indqm;
errstruct(ie).type='dpb2'; 
errstruct(ie).sigma=dkk_quad;
ie=ie+1;

% Sextupoles
%%indsext=[findcells(r0,'Class','Sextupoles')];
indsext=[find(atgetcells(r0,'FamName','SX1')); find(atgetcells(r0,'FamName','SX2')); find(atgetcells(r0,'FamName','SX3'))];
errstruct(ie).indx=indsext;
errstruct(ie).type='psi'; % roll
errstruct(ie).sigma=tilt_sext;
ie=ie+1;
errstruct(ie).indx=indsext;
errstruct(ie).type='x';
errstruct(ie).sigma=dx_sext;
ie=ie+1;
errstruct(ie).indx=indsext;
errstruct(ie).type='y';
errstruct(ie).sigma=dy_sext;
ie=ie+1;
errstruct(ie).indx=indsext;
errstruct(ie).type='s';
errstruct(ie).sigma=ds_sext;
ie=ie+1;

errstruct(ie).indx=indsext;
errstruct(ie).type='dpb3'; 
errstruct(ie).sigma=dkk_sext;
ie=ie+1;

indm=find(atgetcells(r0,'FamName','BPMx'));
sBPM=findspos(ring,indm);
%% set errors

magindex=arrayfun(@(a)a.indx,errstruct,'un',0);
type=arrayfun(@(a)a.type,errstruct,'un',0);
sigma=arrayfun(@(a)a.sigma,errstruct,'un',0);



    
    %% Error loop
    
global GLOBVAL
GLOBVAL.E0=0.05e9;
GLOBVAL.LatticeFile='test';
    
nerr=100;
da_nturns = 50;
dpp = 0.0;

[lindata0, tunes0, chrom0] = twissring(ring_woErr, 0, 1:length(ring_woErr)+1,'chrom', 1e-8); % to get the tunes
beta0=cat(1,lindata0.beta);
SPos=cat(1,lindata0.SPos);
[xxda0,zzda0]=atdynap(ring_woErr,da_nturns,dpp,0.02);
sizebeta=size(beta0);

%
dbeta=zeros([sizebeta nerr]);
Emittance=zeros(nerr,1);
EnergySpread=zeros(nerr,1);
EmitK=zeros(nerr,1);
tunes=zeros(nerr,2);
chrom=zeros(nerr,2);
xorbitmax=zeros(nerr,1);
yorbitmax=zeros(nerr,1);
xorbitrms=zeros(nerr,1);
yorbitrms=zeros(nerr,1);
xorbitmaxErr=zeros(nerr,1);
yorbitmaxErr=zeros(nerr,1);
xorbitrmsErr=zeros(nerr,1);
yorbitrmsErr=zeros(nerr,1);
unstableM = 0;

tic
for kerr=1:nerr
    % Set errors
    fprintf('Trial number          : %8d \n', kerr)

    rerr=atsetrandomerrors(...
    r0,...
    magindex,...
    indm,...
    randi([1 1000],1,1),...
    sigma,...
    2,...
    type);

% define bpm offset and rotation errors
nsigma = 2;
sigma_ox = 200e-6; % random offset errors of 200um
ox=TruncatedGaussian(sigma_ox,nsigma*sigma_ox,[length(indm) 1]);
sigma_oy = 200e-6;
oy=TruncatedGaussian(sigma_oy,nsigma*sigma_oy,[length(indm) 1]);
sigma_gx = 1e-2;% random gain errors of 1-2%
gx=TruncatedGaussian(sigma_gx,nsigma*sigma_gx,[length(indm) 1]);
sigma_gy = 1e-2;
gy=TruncatedGaussian(sigma_gy,nsigma*sigma_gy,[length(indm) 1]);
rx=100e-6; % reading error sigma of 100um (can also be a vector)
ry=100e-6; 
sigma_rot = 1e-5;  % random rotation errors of 10urad
rot=TruncatedGaussian(sigma_rot,nsigma*sigma_rot,[length(indm) 1]);

%dox=1e-5*randn(size(indm)); % random misalignment errors at BPM of 10um
%doy=1e-5*randn(size(indm)); % random misalignment errors at BPM of 10um
% ox=1e-5*randn(size(indm)); % random offset errors of 10um
% oy=1e-5*randn(size(indm)); 
% gx=1e-3*randn(size(indm)); % random gain errors of 0.1%
% gy=1e-3*randn(size(indm));  
% rx=100e-6; % reading error sigma of 100um (can also be a vector)
% ry=100e-6; 
% rot=1e-5*randn(size(indm)); % random rotation errors of 10urad

% set BPM errors
%rerr=atsetshift(rerr,indm,dox,doy);
rerr=atsetbpmerr(rerr,indm,ox,oy,gx,gy,rx,ry,rot);
    
    rerr_struct(:,kerr)  = rerr;
    
    % Get linear data
     
    [lindata, nu, ch] = twissring(rerr, 0, 1:length(rerr)+1,'chrom', 1e-8); % to get the tunes
    beta  =cat(1,lindata.beta);
    alpha = cat(1,lindata.alpha);
    disp  = cat(2,lindata.Dispersion);
    tunes(kerr,:)=nu;
    chrom(kerr,:)=ch;
    %
    dbeta(:,:,kerr)=(beta-beta0)./beta0;
    [Emit, ES, Jx] = atemittance(rerr,beta, alpha, disp);
    Emittance(kerr) = Emit;
    EnergySpread(kerr) = ES;
    EmitK(kerr) = Emit/(1+1/Jx);
    
    % DA
    try
        [xx,zz]=atdynap(rerr,da_nturns,dpp,0.02);
    catch
        warning('Problem using function.  Assigning [xx,zz] value of 0.');
        xx = 0;
        zz = 0;
        unstableM = unstableM + 1;
    end
    xxda(:,kerr)=xx;
    zzda(:,kerr)=zz;
    
    %Orbit
    %orbit4 = findorbit4(rerr,dpp, 1:length(rerr)+1);
    orbit4 = findorbit4(rerr,dpp, indm);
    orbit4Err = findorbit4Err(rerr,dpp, indm);
    
    xorbit(:,kerr)=orbit4(1,:);
    yorbit(:,kerr)=orbit4(3,:);
    xorbitErr(:,kerr)=orbit4Err(1,:);
    yorbitErr(:,kerr)=orbit4Err(3,:);
    
    xorbitmax(kerr,:)=max(abs(orbit4(1,:)));
    yorbitmax(kerr,:)=max(abs(orbit4(3,:)));
    xorbitrms(kerr,:)=rms(orbit4(1,:));
    yorbitrms(kerr,:)=rms(orbit4(3,:));
    
    xorbitmaxErr(kerr,:)=max(abs(orbit4Err(1,:)));
    yorbitmaxErr(kerr,:)=max(abs(orbit4Err(3,:)));
    xorbitrmsErr(kerr,:)=rms(orbit4Err(1,:));
    yorbitrmsErr(kerr,:)=rms(orbit4Err(3,:));
    
end
toc

return
%%
tunesO = tunes;
chromO = chrom;
dbetaO= dbeta;
xxdaO=xxda;
zzdaO=zzda;

xorbitO=xorbit;
yorbitO=yorbit;
xorbitErrO=xorbitErr;
yorbitErrO=yorbitErr;

xorbitmaxO=xorbitmax;
yorbitmaxO=yorbitmax;
xorbitrmsO=xorbitrms;
yorbitrmsO=yorbitrms;

xorbitmaxErrO=xorbitmaxErr;
yorbitmaxErrO=yorbitmaxErr;
xorbitrmsErrO=xorbitrmsErr;
yorbitrmsErrO=yorbitrmsErr;

columnsWithAllZeros = all(xxda == 0);

tunes = tunes(~columnsWithAllZeros,:);
chrom = chrom(~columnsWithAllZeros,:);
dbeta= dbeta(:,:, ~columnsWithAllZeros);
xxda=xxda(:, ~columnsWithAllZeros);
zzda=zzda(:, ~columnsWithAllZeros);
%sum(xxda(:,147)==0)
nonzeroelem = sum(xxda~=0,1);
problem_ind = find(nonzeroelem~=33);
if ~isempty(problem_ind)
% xxda=[xxda(:, 1:problem_ind - 1) xxda(:, problem_ind + 1:end)];
% zzda=[zzda(:, 1:problem_ind - 1) zzda(:, problem_ind + 1:end)];
  index = true(1, size(xxda, 2));
  index([problem_ind]) = false;
  xxda = xxda(:, index);
  zzda = zzda(:, index);
end
xorbit=xorbit(:, ~columnsWithAllZeros);
yorbit=yorbit(:, ~columnsWithAllZeros);
xorbitErr=xorbitErr(:, ~columnsWithAllZeros);
yorbitErr=yorbitErr(:, ~columnsWithAllZeros);

xorbitmax=xorbitmax(~columnsWithAllZeros,:);
yorbitmax=yorbitmax(~columnsWithAllZeros,:);
xorbitrms=xorbitrms(~columnsWithAllZeros,:);
yorbitrms=yorbitrms(~columnsWithAllZeros,:);

xorbitmaxErr=xorbitmaxErr(~columnsWithAllZeros,:);
yorbitmaxErr=yorbitmaxErr(~columnsWithAllZeros,:);
xorbitrmsErr=xorbitrmsErr(~columnsWithAllZeros,:);
yorbitrmsErr=yorbitrmsErr(~columnsWithAllZeros,:);

%%

figure(11);
plot(1e3*xorbitmax,1e3*yorbitmax,'+b','MarkerSize',9,'DisplayName','Without BPM errors')
hold on
plot(1e3*xorbitmaxErr,1e3*yorbitmaxErr,'+r','MarkerSize',9,'DisplayName','With BPM errors')
hold off
set(gcf,'color','w')
set(gca,'fontsize',16');
% xlim([0 8])
% ylim([0 6])
xlim([0 17])
ylim([0 13])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x_{max} [mm]');
ylabel('y_{max} [mm]');
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLmaxErrors_orbitmax_dpp_0015_multip.png','-dpng','-r300')


figure(12);
plot(1e3*xorbitrms,1e3*yorbitrms,'xb','MarkerSize',8,'DisplayName','Without BPM errors')
hold on
plot(1e3*xorbitrmsErr,1e3*yorbitrmsErr,'+r','MarkerSize',8,'DisplayName','With BPM errors')
hold off
set(gcf,'color','w')
set(gca,'fontsize',16');
% xlim([0 4])
% ylim([0 4])
xlim([0 12])
ylim([0 9])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x_{rms} [mm]');
ylabel('y_{rms} [mm]');
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLmaxErrors_orbitrms_dpp_0015_multip.png','-dpng','-r300')


figure(121);
plot(sBPM,1e3.*xorbit(:, 1),'b.-','MarkerSize',10,'DisplayName','Orbit')
hold on; 
plot(sBPM,1e3.*xorbitErr(:, 1),'rx-','MarkerSize',10,'DisplayName','BPM readings');
plot(sBPM,1e3.*xorbit(:, 1:3),'b.-','MarkerSize',10,'HandleVisibility','off')
plot(sBPM,1e3.*xorbitErr(:, 1:3),'rx-','MarkerSize',10,'HandleVisibility','off');
hold off
set(gcf,'color','w')
set(gca,'fontsize',16');
 xlim([0 18])
% ylim([0 6])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('s-position [m]');
ylabel('x [mm]');
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLmaxErrors_xorbit_dpp_0_multip.png','-dpng','-r300')

figure(122);
plot(sBPM,1e3.*yorbit(:, 1),'b.-','MarkerSize',10,'DisplayName','Orbit')
hold on; 
plot(sBPM,1e3.*yorbitErr(:, 1),'rx-','MarkerSize',10,'DisplayName','BPM readings');
plot(sBPM,1e3.*yorbit(:, 1:3),'b.-','MarkerSize',10,'HandleVisibility','off')
plot(sBPM,1e3.*yorbitErr(:, 1:3),'rx-','MarkerSize',10,'HandleVisibility','off');
hold off
set(gcf,'color','w')
set(gca,'fontsize',16');
 xlim([0 18])
% ylim([0 6])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('s-position [m]');
ylabel('y [mm]');
u = legend('show','Location','NorthEast');
set(u, 'FontSize',14)
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLmaxErrors_yorbit_dpp_0_multip.png','-dpng','-r300')



dbetax=squeeze(dbeta(:,1,:));sdbetax=std(dbetax');maxdbetax=max(dbetax');dbetaxmax=max(abs(dbetax));
dbetaz=squeeze(dbeta(:,2,:));sdbetaz=std(dbetaz');maxdbetaz=max(dbetaz');dbetazmax=max(abs(dbetaz));

figure(13);
plot(SPos,sdbetax*100,'-r', 'Linewidth', 1.6);hold on
plot(SPos,sdbetaz*100,'-b', 'Linewidth', 1.6);hold off
axis tight
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
legend('Horizontal','Vertical')
xlabel('s-position [m]');               
ylabel('[\Delta\beta/\beta]_{rms} [%]')
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLmaxErrors_betabeatrms_vs_s_multip.png','-dpng','-r300')

figure(14);
plot(SPos,maxdbetax*100,'-r', 'Linewidth', 1.6);hold on
plot(SPos,maxdbetaz*100,'-b', 'Linewidth', 1.6);hold off
axis tight
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
legend('Horizontal','Vertical')
xlabel('s-position [m]');               
ylabel('[\Delta\beta/\beta]_{max} [%]')
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLmaxErrors_betabeatmax_vs_s_multip.png','-dpng','-r300')

figure(15);
plot(dbetaxmax,dbetazmax,'ro');
axis tight
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('[\Delta\beta_x/\beta_x]_{max}');
ylabel('[\Delta\beta_z/\beta_z]_{max}');
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLmaxErrors_betabeatmax_dpp_0_multip.png','-dpng','-r300')

figure(16);
plot(tunes(:,1),tunes(:,2),'ob','MarkerSize',8);
hold on
plot(tunes0(:,1),tunes0(:,2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',9)
hold off
axis tight
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
u=legend('Random seeds','Initial tunes (3.17/1.64)')
set(u, 'Location','NorthEast')
xlabel('\nu_x');
ylabel('\nu_z');
addlabel(1, 0, datestr(clock,0))
%print('thomx_ALLmaxErrors_tunes_dpp_0_multip.png','-dpng','-r300')


%% Save data

E.file='ThomX_017_064_r56_02_chro00_AT2';

E.nseeds = nerr;
E.ringerrors = rerr_struct;

E.dipole_x=dx_dipole;
E.dipole_y=dy_dipole;
E.dipole_s=ds_dipole;
E.dipole_tilt=tilt_dipole;
E.dkk_dipole=dkk_dipole;

E.quad_x=dx_quad;
E.quad_y=dy_quad;
E.quad_s=ds_quad;
E.quad_tilt=tilt_quad;
E.dkk_quad=dkk_quad;

E.sext_x=dx_sext;
E.sext_y=dy_sext;
E.sext_s=ds_sext;
E.sext_tilt=tilt_sext;
E.dkk_sext=dkk_sext;

E.BPM_oxy=ox;
E.BPM_gxy=gx;
E.BPM_rxy=rx;
E.BPM_rot=rot;

E.at=r0;
E.SPos=SPos;
E.nerr=nerr;
E.beta0=beta0;
E.tunes0=tunes0;
E.dbeta=dbeta;
E.tunes=tunes;
E.emittance=Emittance;

E.da_nturns = da_nturns;
E.xda0=xxda0;
E.zda0=zzda0;
E.xda=xxda;
E.zda=zzda;

E.xorbitmax = xorbitmax;
E.yorbitmax = yorbitmax;
E.xorbitrms = xorbitrms;
E.yorbitrms = yorbitrms;

E.xorbitmaxErr = xorbitmaxErr;
E.yorbitmaxErr = yorbitmaxErr;
E.xorbitrmsErr = xorbitrmsErr;
E.yorbitrmsErr = yorbitrmsErr;

E.xorbit = xorbit;
E.yorbit= yorbit;
E.xorbitErr = xorbitErr;
E.yorbitErr = yorbitErr;

%errors max
%save data_dipole_align_err_100um_05mrad E errstruct
%save data_dipole_field_err_0005 E errstruct
%save data_quad_align_err_100um_05mrad E errstruct
%save data_quad_field_err_0005 E errstruct
%save data_sext_align_err_100um_05mrad E errstruct
%save data_sext_field_err_0005 E errstruct
% save data_ALLerrorsDA_align_err_100um_05mrad_field_err_0005 E errstruct

%errors min
%save data_dipole_align_err_30um_02mrad E errstruct
%save data_dipole_field_err_0001 E errstruct
%save data_quad_align_err_30um_02mrad E errstruct
%save data_quad_field_err_0001 E errstruct
%save data_sext_align_err_30um_02mrad E errstruct
%save data_sext_field_err_0001 E errstruct
%save data_ALLerrorsDA_align_err_30um_02mrad_field_err_0001 E errstruct

%DA off momentum
% save data_ALLerrorsDA_align_err_30um_02mrad_field_err_0001_dpp_001 E errstruct
%save data_ALLerrorsDA_align_err_30um_02mrad_field_err_0001_dpp_0 E errstruct

%Multipole+BPM error off 
%save data_ALLerrorsMultipBPM_align_err_30um_02mrad_field_err_0001_dpp_0 E errstruct
%save data_ALLerrorsMultipBPM_align_err_100um_05mrad_field_err_0005_QUAD30um0mrad0_dpp_0 E errstruct
%save data_ALLerrorsMultipBPM_align_err_100um_05mrad_field_err_0005_dpp_0_GOOD E errstruct

%DA off momentum
%save data_ALLerrorsMultipBPM_align_err_30um_02mrad_field_err_0001_dpp_0015 E errstruct
%save data_ALLerrorsMultipBPM_align_err_30um_02mrad_field_err_0001_dpp_m0015 E errstruct
%save data_ALLerrorsMultipBPM_align_err_30um_02mrad_field_err_0001_dpp_0_100seed E errstruct

%save data_ALLerrorsMultipBPM_align_err_100um_05mrad_field_err_0005_dpp_m001_GOOD E errstruct
save data_ALLerrorsMultipBPM_align_err_100um_05mrad_field_err_0005_dpp_0_GOOD_100seed E errstruct


%% DA


figure(17);
plot(xxda0*1e3,zzda0*1e3,'b-','LineWidth',2);hold on,
plot(xxda*1e3,zzda*1e3,'color',[0,0,0]+0.5,'LineWidth',1);
plot(xxda0*1e3,zzda0*1e3,'b-','LineWidth',2);hold off
legend('Initial: no errors','Trials with errors')
% xlim([-50 50])
% ylim([0 30])
xlim([-30 30])
ylim([0 20])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x [mm]');                 % Add labels
ylabel('z [mm]');
title('DA (\deltap = 0%)')
print('thomx_ALLmaxErrors_DA_init_seeds_dpp_0_multip100.png','-dpng','-r300')

% Stat over DA
xxdai=[-40:0.2:40]*1e-3;
zzdai0=interp1(xxda0,zzda0,xxdai','pchip',0);
for kerr=1:length(xxda)%nerr
    zzdai(:,kerr)=interp1(xxda(:,kerr),zzda(:,kerr),xxdai','pchip',0);
end

mzzdai=mean(zzdai,2);
szzdai=std(zzdai')';
DA_surf=sum(zzdai)*0.1e-3; % surface in m^2
DA_surf0=sum(zzdai0)*0.1e-3; % surface in m^2

figure(18)
plot(xxdai*1e3,zzdai0*1e3,'-b','LineWidth',2);hold on
plot(xxdai*1e3,mzzdai*1e3,'-r','LineWidth',2);
plot(xxdai*1e3,(mzzdai-szzdai)*1e3,'--m','LineWidth',2);
plot(xxdai*1e3,(mzzdai+szzdai)*1e3,'--m','LineWidth',2);hold off
legend('Initial: no errors','Mean of 100 seeds','Mean - \sigma','Mean + \sigma')
% xlim([-60 60])
% ylim([0 25])
xlim([-30 30])
ylim([0 20])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x [mm]');                 % Add labels
ylabel('z [mm]');
title('DA (\deltap = 0%)')
print('thomx_ALLmaxErrors_DA_init_,meanseeds1_dpp_0_multip100.png','-dpng','-r300')

figure(181)
plot(xxda0*1e3,zzda0*1e3,'b-','LineWidth',2);hold on
plot(xxdai*1e3,mzzdai*1e3,'-r','LineWidth',2);
plot(xxda*1e3,zzda*1e3,'color',[0,0,0]+0.5,'LineWidth',1);
plot(xxdai*1e3,mzzdai*1e3,'-r','LineWidth',2);
plot(xxda0*1e3,zzda0*1e3,'b-','LineWidth',2)
% plot(xxda0*1e3,zzda0*1e3,'m-','LineWidth',2)
%plot(xxdai*1e3,(mzzdai-szzdai)*1e3,'--m','LineWidth',2);
%plot(xxdai*1e3,(mzzdai+szzdai)*1e3,'--m','LineWidth',2);hold off
legend('Initial: no errors','Mean of 100 seeds')
xlim([-30 30])
ylim([0 17])
grid on
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('x [mm]');                 % Add labels
ylabel('z [mm]');
title('DA (\deltap = 0%)')
print('thomx_ALLmaxErrors_DA_init_meanseeds2_dpp_0_multip100.png','-dpng','-r300')

figure(19);
plot(DA_surf0*1e6, 0 ,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10);hold on
h=histogram(DA_surf*1e6); 
h.FaceColor = [0 0.5 0.5];
hold off    
set(gcf,'color','w')
set(gca,'fontsize',16');
xlabel('DA surface [mm^2]');
ylabel('Entries');
title('DA (\deltap = 0%)')
u=legend('Initial: no errors','Trials with errors')
set(u,'Location','NorthEast','Orientation','vertical')
print('thomx_ALLmaxErrors_DA_surf_dpp_0_multip100.png','-dpng','-r300')


%%

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
atplot(rerr,'comment',[],@plClosedOrbit)


figure('units','normalized','position',[0.1 0.4 0.45 0.35])
atplot(rerr,'comment',[],@pltmisalignments);
%print('thomx_ALLmaxErrors_misalignments_dpp_0_multip.png','-dpng','-r300')
%export_fig('thomx_ALLminErrors_misalignments_dpp_0_multip.pdf','-r300');

return

%% Check the bad machines

bad_machines = find(columnsWithAllZeros~=0);
indq=find(atgetcells(ring,'Class','Quadrupole'));
indsext=[find(atgetcells(r0,'FamName','SX1')); find(atgetcells(r0,'FamName','SX2')); find(atgetcells(r0,'FamName','SX3'))];
inds=findcells(r0,'Class','bend');


for ib = 1:length(bad_machines)
[X,Y,S,Tilt]=GetMisalignments(rerr_struct(:,bad_machines(ib)))
%PolynomBVal1=getcellstruct(ring,'PolynomB',inds,1,1);
PolynomBVal1err=getcellstruct(rerr_struct(:,bad_machines(ib)),'PolynomBErr',inds);

PolynomBVal2=getcellstruct(ring,'PolynomB',indq,1,2);
PolynomBVal2err=getcellstruct(rerr_struct(:,bad_machines(ib)),'PolynomB',indq,1,2);

PolynomBVal3=getcellstruct(ring,'PolynomB',indsext,1,3);
PolynomBVal3err=getcellstruct(rerr_struct(:,bad_machines(ib)),'PolynomB',indsext,1,3);

FieldErr_DIP(:,ib) = PolynomBVal1err;
FieldErr_QUAD(:,ib) = (PolynomBVal2 -PolynomBVal2err)./PolynomBVal2;
FieldErr_SEXT(:,ib) = (PolynomBVal3 -PolynomBVal3err)./PolynomBVal3;

X_DIP(:,ib) = X(inds)
Y_DIP(:,ib) = Y(inds)
S_DIP(:,ib) = S(inds)
Tilt_DIP(:,ib) = Tilt(inds)

X_QUAD(:,ib) = X(indq)
Y_QUAD(:,ib) = Y(indq)
S_QUAD(:,ib) = S(indq)
Tilt_QUAD(:,ib) = Tilt(indq)

X_SEXT(:,ib) = X(indsext)
Y_SEXT(:,ib) = Y(indsext)
S_SEXT(:,ib) = S(indsext)
Tilt_SEXT(:,ib) = Tilt(indsext)
end


for ib = 1:100
%indq=find(atgetcells(ring,'Class','Quadrupole'))
[X,Y,S,Tilt]=GetMisalignments(rerr_struct(:,ib))

%PolynomBVal1=getcellstruct(ring,'PolynomB',inds,1,1);
PolynomBVal1err=getcellstruct(rerr_struct(:,ib),'PolynomBErr',inds);

PolynomBVal2=getcellstruct(ring,'PolynomB',indq,1,2);
PolynomBVal2err=getcellstruct(rerr_struct(:,ib),'PolynomB',indq,1,2);

PolynomBVal3=getcellstruct(ring,'PolynomB',indsext,1,3);
PolynomBVal3err=getcellstruct(rerr_struct(:,ib),'PolynomB',indsext,1,3);

FieldErr_DIP_GM(:,ib) = PolynomBVal1err;
FieldErr_QUAD_GM(:,ib) = (PolynomBVal2 -PolynomBVal2err)./PolynomBVal2;
FieldErr_SEXT_GM(:,ib) = (PolynomBVal3 -PolynomBVal3err)./PolynomBVal3;

X_DIP_GM(:,ib) = X(inds)
Y_DIP_GM(:,ib) = Y(inds)
S_DIP_GM(:,ib) = S(inds)
Tilt_DIP_GM(:,ib) = Tilt(inds)

X_QUAD_GM(:,ib) = X(indq)
Y_QUAD_GM(:,ib) = Y(indq)
S_QUAD_GM(:,ib) = S(indq)
Tilt_QUAD_GM(:,ib) = Tilt(indq)

X_SEXT_GM(:,ib) = X(indsext)
Y_SEXT_GM(:,ib) = Y(indsext)
S_SEXT_GM(:,ib) = S(indsext)
Tilt_SEXT_GM(:,ib) = Tilt(indsext)
end


%%
figure
plot(max(X_QUAD_GM), max(Y_QUAD_GM),'bo','MarkerSize',8)
hold on
plot(max(X_QUAD), max(Y_QUAD),'ro','MarkerSize',8)
hold off
xlabel('X');
ylabel('Y');

figure
plot(max(X_DIP_GM), max(Y_DIP_GM),'bo','MarkerSize',8)
hold on
plot(max(X_DIP), max(Y_DIP),'go','MarkerSize',8)
hold off
xlabel('X');
ylabel('Y');

figure
plot(max(X_SEXT_GM), max(Y_SEXT_GM),'bo','MarkerSize',8)
hold on
plot(max(X_SEXT), max(Y_SEXT),'mo','MarkerSize',8)
hold off
xlabel('X');
ylabel('Y');

%%

figure
plot(max(S_QUAD_GM), max(Tilt_QUAD_GM),'bo','MarkerSize',8)
hold on
plot(max(S_QUAD), max(Tilt_QUAD),'ro','MarkerSize',8)
hold off
xlabel('S');
ylabel('Tilt');

figure
plot(max(Tilt_DIP_GM),'bo','MarkerSize',8)
hold on
plot(max(Tilt_DIP),'go','MarkerSize',8)
hold off
%xlabel('S');
ylabel('Tilt');

figure
plot(max(S_SEXT_GM), max(Tilt_SEXT_GM),'bo','MarkerSize',8)
hold on
plot(max(S_SEXT), max(Tilt_SEXT),'mo','MarkerSize',8)
hold off
xlabel('S');
ylabel('Tilt');

%%

figure
plot(max(FieldErr_QUAD_GM),'bo','MarkerSize',8)
hold on
plot(max(FieldErr_QUAD),'ro','MarkerSize',8)
hold off
ylabel('Field Error');

figure
plot(max(FieldErr_DIP_GM),'bo','MarkerSize',8)
hold on
plot(max(FieldErr_DIP),'go','MarkerSize',8)
hold off
ylabel('Field Error');

figure
plot(max(FieldErr_SEXT_GM),'bo','MarkerSize',8)
hold on
plot(max(FieldErr_SEXT),'mo','MarkerSize',8)
hold off
ylabel('Field Error');


