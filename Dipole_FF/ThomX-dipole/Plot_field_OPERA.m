%function  Plot_field_OPERA
%PLOT_FIELD Summary of this function goes here
%   Detailed explanation goes here
%   Maille ThomX
% Hard hedge traj

%Dipole 
radius=0.352; %m
teta  =pi/8;    % half field deviation 45??
tetapole=pi/9;  % half yoke reduce to 40??
gap   =0.042 ;
fprintf('Some plot from OPERA mTm \n')

%[Xb ,Sb , Bz] = getOPERAfield2D('field/champ_dipole_OPERA.table');
%[Xb ,Sb , Bz] = getOPERAfield2D('field/champ_dipole_OPERA_v6.table');
%[Xb ,Sb , Bz] = getOPERAfield2D('field/champ_dipole_OPERA_v31.table');
% [Xb ,Sb , Bz] = getOPERAfield2D('field/champ_dipole_OPERA_v33bis.table');
%[Xb ,Sb , Bz] = getOPERAfield2D('field/dipole33_XC10_2016_200A.table');
[Xb ,Sb , Bz] = getOPERAfield2D('field/dipole33_XC10.table');

%[Xb,Sb,Bz] = getMEASUREDfield2D('field/20160204c_THOMX#009_fieldmap_200A.xlsx');
%[Xb,Sb,Bz] = getMEASUREDfield2D('field-meas/20160129a_THOMX#009_fieldmap_300A.xlsx');
Bz=-Bz*0.7965; % *0.7965 pour ~50 MeV
bz0=-interp2(Xb,Sb,Bz,0,0);
fprintf('  Peak field  %g  T \n',bz0 )


[pole] = make_pole(radius,0.050,2*tetapole);
figure(1)
set(gcf,'color','w')
set(gca,'fontsize',18)
plot(pole(2,:),pole(1,:),'-w'); 
hold on
image(Sb(:,1),Xb(1,:),-Bz'*100)
plot(pole(2,:),pole(1,:),'-w','Linewidth',3); 
xlim([-250 250])
ylim([-100 100])
xlabel('ds (mm)')
ylabel('dx (mm)')
grid on
hold off
colorbar
print('dipole_img','-dpng','-r300')

figure(11)
set(gcf,'color','w')
set(gca,'fontsize',18)
plot(pole(2,:),pole(1,:),'-w'); 
hold on
imagesc(Sb(:,1),Xb(1,:),-Bz');
plot(pole(2,:),pole(1,:),'-w','Linewidth',3); 
shading interp;
xlim([-250 250])
ylim([-100 100])
xlabel('s (mm)')
ylabel('x (mm)')
grid on
hold off
colorbar
print('dipole_imgsc','-dpng','-r300')

   
   
% 

%     fprintf(' \n')
%     fprintf('  Fit at center: \n')
%     
%     % Simple Bz on the middle
%     xx =Xb(1,81:121);xx1=xx;
%     dbz=(-Bz(500,81:121)-bz0)/bz0;
%     p = polyfit(xx*1e-3,dbz,4);
%     fprintf('    %8i',[4:-1:0]);fprintf('\n')
%     fprintf('    %8.2e',p);fprintf('\n')
%     fprintf('    %8.2e',p(1)*(2e-2)^4,p(2)*(2e-2)^3,...
%                         p(3)*(2e-2)^2,p(4)*(2e-2)^1);fprintf('\n')
%     f = polyval(p,xx*1e-3); 
%     figure(5)
%     set(gca,'Fontsize',14)
%     plot(xx,dbz*100,'ob'); hold on
%     plot(xx,f*100,'--b'); hold off
%     xlim([-21 21])
%     xlabel('dx (mm)')
%     ylabel('dB/B (%)')
%     legend('Radial field profil', 'Fit')
%     grid on
%     fprintf(' \n')

%
[ traj0 , ds ]      = make_traj(radius,0,teta);
[ trajp20 , dsp20 ] = make_traj((radius+0.020),+0.020,teta);
[ trajm20 , dsm20 ] = make_traj((radius-0.020),-0.020,teta);
[pole] = make_pole(radius,0.050,2*tetapole);


%
figure(6)
set(gcf,'color','w')
set(gca,'fontsize',18)
plot(traj0(2,:),traj0(1,:),'-k');hold on % trick to get image in good direction
%image(Sb(:,1),Xb(1,:),-Bz'*100); 
imagesc(Sb(:,1),Xb(1,:),-Bz');
plot(pole(2,:),pole(1,:),'-w','Linewidth',3); 
plot(traj0(2,:),traj0(1,:),'-w','Linewidth',2);
plot(trajp20(2,:),trajp20(1,:),'-w','Linewidth',2);
plot(trajm20(2,:),trajm20(1,:),'-w','Linewidth',2);hold off
shading interp;
colorbar
xlim([0 300]);
ylim([-100 30 ])
xlabel('s (mm)')
ylabel('x (mm)')
grid on
print('dipole_traj','-dpng','-r300')


%%


%%

%
bz   =interp2(Xb,Sb,Bz,traj0(1,:)  ,traj0(2,:));
bzp20=interp2(Xb,Sb,Bz,trajp20(1,:),trajp20(2,:));
bzm20=interp2(Xb,Sb,Bz,trajm20(1,:),trajm20(2,:));
intbz0=sum(-bz.*ds)*1000;
fprintf('  Integrated field  %g  mTm \n',intbz0 *2)
fprintf('  Equivalent length %g  mm \n',-intbz0/bz(1)*2)
K= sum(-bz.*(bz0+bz).*ds)/bz0^2/gap;
fprintf('  Fringe coef K = %g  \n',K)

figure(7)
set(gca,'Fontsize',18)
plot(cumsum(ds)*1000   ,-bz,'-k');hold on
plot(cumsum(dsp20)*1000,-bzp20,'-b');
plot(cumsum(dsm20)*1000,-bzm20,'-r');hold off
legend('dx=0','dx=+20mm','dx=-20mm')
xlabel('s-traj (mm)')
ylabel('B (T)')
grid on
print('dipole_fieldtraj','-dpng','-r300')

% s=cumsum(ds);
% save ThomX_dipole_field_50MeV.mat s bz

% integrated profil
xx=(-20:2:20);
intbz=[];
for dx=xx
    [ traj , ds ] = make_traj((radius+dx*1e-3),dx*1e-3,teta);
    bz=interp2(Xb,Sb,Bz,traj(1,:),traj(2,:));
    intbz=[intbz sum(-bz.*ds)*1000];
end

fprintf(' \n')
fprintf('  Integrated: \n')
intbzr =(intbz-intbz0)/intbz0;
intbzrc=(intbz-intbz0 - bz0*tetapole*xx)/intbz0;
p = polyfit(xx*1e-3,intbzrc,4);pint=p;
%p1 =[p(1)*0 p(2)*0 p(3) p(4)*0];
fprintf('    %8i',[4:-1:0]);fprintf('\n')
fprintf('    %8.2e',p);fprintf('\n')
fprintf('    %8.2e',p(1)*(2e-2)^4,p(2)*(2e-2)^3,...
                    p(3)*(2e-2)^2,p(4)*(2e-2)^1);fprintf('\n')
f = polyval(p,xx*1e-3); 
figure(9)
set(gca,'Fontsize',14)
plot(xx,intbzr*100,'ob');hold on
plot(xx,intbzrc*100,'or');
plot(xx,f*100,'--r'); hold off
xlabel('dx (mm)')
ylabel('d intB / intB (%)')
legend('Direct int field profil','Corr int field profil', 'Fit')
grid on

% figure(10)
% set(gca,'Fontsize',14)
% plot(xx,intbzrc*100,'or');hold on
% plot(xx,f*100,'--r'); hold off
% xlabel('dx (mm)')
% ylabel('d intB / intB (%)')
% legend('Corr int field profil', 'Fit')
% grid on


fprintf(' \n')







