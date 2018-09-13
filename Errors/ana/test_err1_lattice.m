%% test lattice with error fields
clear
global GLOBVAL
GLOBVAL.E0=0.05e9;
GLOBVAL.LatticeFile='test';
%dir ='/home/sources/physmach/loulergue/work/matlab/Simulation/beta_structure/';
% 125 pm
% input_file =[dir 'Upgrade/test-7BA/centre_opt1_7BA_type1_fit9_atmatch3-refit1.str'];
% RING=atreadbeta_alex(input_file);
% RING=fitbend(RING,pi/8,[1.02 1.0 1.0 1.0]*1, '7BAsym16');
% RING=fittune2_alex(RING, [2.763 0.888], 'QF68E', 'B61HE');
% RING=fitchrom_alex(RING,[-0.1 -0.12],'SXD1E' ,'SXF1E' );
% RING=addoct(RING,[36  -69  30],'7BA');

% 72 pm 
%load Lattice-save/SOLEIL_U_v9 

RING0 = ThomX_017_064_r56_02_chro00;
% Develop over cells
% RING0=[];
% for k=1:20
%     RING0=[RING0 ; RING];
% end

%[lindata0, tunes0, chrom0]=atlinopt(RING0,0,1:length(RING0)+1);
[lindata0, tunes0, chrom0] = twissring(RING0, 0, 1:length(RING0)+1,'chrom', 1e-8); % to get the tunes
beta0=cat(1,lindata0.beta);
SPos=cat(1,lindata0.SPos);
[xxda0,zzda0]=atdynap(RING0,25,0.0,0.02);
sizebeta=size(beta0);

%% Error loop
nerr=200;
dksk=1e-3;
%
dbeta=zeros([sizebeta nerr]);
Emittance=zeros(nerr,1);
EnergySpread=zeros(nerr,1);
EmitK=zeros(nerr,1);
tunes=zeros(nerr,2);
chrom=zeros(nerr,2);
for kerr=1:nerr
    % Set errors
    fprintf('Trial number          : %8d \n', kerr)
    [RING]=set_B2_error(RING0,'Bend',dksk);
    [RING]=set_B2_error(RING,'Quadrupole',dksk);
    % Get linear data
    %[lindata, tunes, chrom]=atlinopt(RING,0,1:length(RING)+1);
    [lindata, nu, ch] = twissring(RING, 0, 1:length(RING)+1,'chrom', 1e-8); % to get the tunes
    beta  =cat(1,lindata.beta);
    alpha = cat(1,lindata.alpha);
    disp  = cat(2,lindata.Dispersion);
    tunes(kerr,:)=nu;
    chrom(kerr,:)=ch;
    %
    dbeta(:,:,kerr)=(beta-beta0)./beta0;
    [Emit, ES, Jx] = atemittance(RING,beta, alpha, disp);
    Emittance(kerr) = Emit;
    EnergySpread(kerr) = ES;
    EmitK(kerr) = Emit/(1+1/Jx);
    % DA
    [xx,zz]=atdynap(RING,25,0,0.02);
    xxda(:,kerr)=xx;
    zzda(:,kerr)=zz;
end


% plot
figure(1);
dbetax=squeeze(dbeta(:,1,:));sdbetax=std(dbetax');
dbetaz=squeeze(dbeta(:,2,:));sdbetaz=std(dbetaz');
set(gcf,'color','w')
set(gca,'fontsize',16');
plot(SPos,sdbetax*100,'-r');hold on
plot(SPos,sdbetaz*100,'-b');hold off
axis tight
legend('H','V')
xlabel('S [m]');               
ylabel('rms beta beat [%]');

figure(11);
set(gcf,'color','w')
set(gca,'fontsize',16');
plot(tunes0(:,1),tunes0(:,2),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10);hold on
plot(tunes(:,1),tunes(:,2),'ob');hold off
axis tight
legend('Initial','Trials')
xlabel('\nu_x');
ylabel('\nu_z');

figure(2);
set(gcf,'color','w')
set(gca,'fontsize',16');
hist(Emittance*1e12);             
xlabel('Emittance [pm]');

figure(3);
set(gcf,'color','w')
set(gca,'fontsize',16');
plot(xxda0*1e3,zzda0*1e3,'-k','LineWidth',2);hold on,
plot(xxda*1e3,zzda*1e3,'-b','LineWidth',1);
plot(xxda0*1e3,zzda0*1e3,'-k','LineWidth',2);hold off
legend('Initial','Trials')
xlim([-50 50])
ylim([0 25])
grid on
xlabel('X [mm]');                 % Add labels
ylabel('Z [mm]');
title('DA (\delta = 0)')

% Stat over DA
xxdai=[-50:0.1:50]*1e-3;
zzdai0=interp1(xxda0,zzda0,xxdai','cubic',0);
for kerr=1:nerr
    zzdai(:,kerr)=interp1(xxda(:,kerr),zzda(:,kerr),xxdai','cubic',0);
end

mzzdai=mean(zzdai,2);
szzdai=std(zzdai')';
DA_surf=sum(zzdai)*0.1e-3; % surface in m^2
DA_surf0=sum(zzdai0)*0.1e-3; % surface in m^2

figure(31)
set(gcf,'color','w')
set(gca,'fontsize',16');
plot(xxdai*1e3,zzdai0*1e3,'-k','LineWidth',2);hold on
plot(xxdai*1e3,mzzdai*1e3,'-b','LineWidth',2);
plot(xxdai*1e3,(mzzdai-szzdai)*1e3,'-g','LineWidth',2);
plot(xxdai*1e3,(mzzdai+szzdai)*1e3,'-g','LineWidth',2);hold off
legend('Initial','Mean','Mean - 1*sig')
xlim([-50 50])
ylim([0 25])
grid on
xlabel('X [mm]');                 % Add labels
ylabel('Z [mm]');
title('DA (\delta = 0)')

figure(4);
set(gcf,'color','w')
set(gca,'fontsize',16');
plot(DA_surf0*1e6, 0 ,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10);hold on
hist(DA_surf*1e6); hold off            
xlabel('DA surface [mm^2]');
legend('Initial','Trials')

return
%%
E.file='ThomX_017_064_r56_02_chro00';
E.at=RING0;
E.SPos=SPos;
E.bend=dksk;
E.quad=dksk;
E.nerr=nerr;
E.beta0=beta0;
E.tunes0=tunes0;
E.xda0=xxda0;
E.zda0=zzda0;
E.dbeta=dbeta;
E.tunes=tunes;
E.emittance=Emittance;
E.xda=xxda;
E.zda=zzda;

%file='data_err1_125pm_B2_1permil'
save data_err1_75pm_B2_1permil E

%%
[xda0,zda0]=atdynap(RING0,25,0.0,0.02);

return
