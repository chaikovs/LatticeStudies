%% test lattice with error fields from saved data by test_err_lattice
clear
%load data_err1_125pm_B2_1permil E
load data_err1_75pm_B2_1permil E


% base data
RING0=E.at;
SPos=E.SPos;
dksk=E.bend;
nerr=E.nerr;
beta0=E.beta0;
tunes0=E.tunes0;
xxda0=E.xda0;
zzda0=E.zda0;
% data over nerr trials
dbeta=E.dbeta;
tunes=E.tunes;
Emittance=E.emittance;
xxda=E.xda;
zzda=E.zda;

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
xlim([-5 5]*0.9)
ylim([0 5])
grid on
xlabel('X [mm]');                 % Add labels
ylabel('Z [mm]');
title('DA (\delta = 0)')

% Stat over DA
xxdai=[-5:0.1:5]*1e-3;
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
xlim([-5 5]*0.9)
ylim([0 5])
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


