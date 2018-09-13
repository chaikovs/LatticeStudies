%function  Plot_field_OPERA
%PLOT_FIELD Summary of this function goes here
%   Detailed explanation goes here
%   Maille ThomX
% Hard hedge traj

%Dipole 
clear
radius  =0.352; %m
teta    =pi/8; % half field deviation
tetapole=pi/9; % half yoke deviation
gap     =0.042 ;

fprintf(' \n')
fprintf('######################## \n')
fprintf('Some plot from OPERA mTm \n')

%[Xb ,Sb , Bz] = getOPERAfield2D('field/champ_dipole_OPERA.table');
%[Xb ,Sb , Bz] = getOPERAfield2D('field/champ_dipole_OPERA_v6.table');
%[Xb ,Sb , Bz] = getOPERAfield2D('field/champ_dipole_OPERA_v15.table');
%[Xb ,Sb , Bz] = getOPERAfield2D('field/champ_dipole_OPERA_v33_CC.table');
%[Xb ,Sb , Bz] = getOPERAfield2D('field/dipole-perfect_XC10.table');
%[Xb ,Sb , Bz] = getOPERAfield2D('field/dipole-defaut050??m_XC10.table');
[Xb,Sb,Bz] = getMEASUREDfield2D('field-meas/20160129a_THOMX#009_fieldmap_300A.xlsx');
bz0=interp2(Xb,Sb,Bz,0,0);
fprintf('   Peak field  %g  T \n',bz0 )

%
[ traj0 , ds , np ] = make_traj(radius,0,teta);
bz   =interp2(Xb,Sb,Bz,traj0(1,:),traj0(2,:));
ds0=ds(1);
Lmag  =sum(bz.*ds)/bz(1); % equivalent half length

bz   =interp2(Xb,Sb,Bz,traj0(1,:)  ,traj0(2,:));
% bzp20=interp2(Xb,Sb,Bz,trajp20(1,:),trajp20(2,:));
% bzm20=interp2(Xb,Sb,Bz,trajm20(1,:),(2,:));
intbz0=sum(bz.*ds);

% Get profil along path bzn
xx=(-20:2:20);
intbz=[];
bzn=[];
for dx=xx
    [ traj , ds ] = make_traj((radius+dx*1e-3),dx*1e-3,teta);
    bz=interp2(Xb,Sb,Bz,traj(1,:),traj(2,:));
    intbz=[intbz sum(bz.*ds)];
    bzn=[bzn ; bz.*ds/ds0];
    %bzn=[bzn ; bz];
end


% Get poly fit on integrated form
intbzr =(intbz-intbz0)/intbz0;
intbzrc=(intbz-intbz0 - bz0*tetapole*xx/1000)/intbz0;
p = polyfit(xx*1e-3,intbzrc,4);pint=p;

% Get polyfit along path 
pn=[];
for i=1:1:length(ds)
    corr=0;% to remove natural dip H foc in matrix
    if i<(np+2) ; corr=bzn(11,i)*xx'*1e-3*teta/(np+1)/ds0;end
    db=(bzn(:,i)-bzn(11,i)-corr)/bzn(11,1);
    p = polyfit(xx*1e-3,db',4);
    pn=[pn ; p];
end

% get polyfit sum
% sum bit diff from above : ds variation not included !
ps=sum(pn)*ds0/Lmag;

%load the comp dipole v10
%load dipole_v10.mat traj0v10 pnv10
% load dipole_v16.mat traj0v16 pnv16
% traj0v10= traj0v16;pnv10=pnv16;
%  load dipole_v19.mat traj0v19 pnv19
%  traj0v10= traj0v19;pnv10=pnv19;
%  load dipole_v21.mat traj0v21 pnv21
%  traj0v10= traj0v21;pnv10=pnv21;
load dipole_v33bis.mat traj0v33bis pnv33bis
traj0v10= traj0v33bis;pnv10=pnv33bis;
% load dipole_perfect.mat traj0perfect pnvperfect
% traj0v10= traj0perfect;pnv10=pnvperfect;
 
% plot grad
figure(2)
set(gcf,'color','w')
set(gca,'fontsize',16)
plot(traj0v10(2,:),pnv10(:,4),'-b','linewidth',2);hold on
plot(traj0(2,:),pn(:,4),'-r','linewidth',2);hold off
xlabel('S (mm)')
ylabel('Gradient')
legend('Simu OPERA','Meas @ALBA','Location','Southeast')
grid on
fprintf('   Sum quad= %8.2e     int %8.2e\n',ps(4), pint(4))

% plot sext
figure(3)
set(gcf,'color','w')
set(gca,'fontsize',16)
plot(traj0v10(2,:),pnv10(:,3),'-b','linewidth',2);hold on
plot(traj0(2,:),pn(:,3),'-r','linewidth',2);hold off
xlabel('S (mm)')
ylabel('Sextupole')
legend('Simu OPERA','Meas ALBA','Location','Southeast')
grid on
fprintf('   Sum sext= %8.2e     int %8.2e \n',ps(3),pint(3))

% plot oct
figure(4)
set(gcf,'color','w')
set(gca,'fontsize',16)
plot(traj0v10(2,:),pnv10(:,2),'-b','linewidth',2);hold on
plot(traj0(2,:),pn(:,2),'-r','linewidth',2);hold off
xlabel('S (mm)')
ylabel('Octupole')
legend('Simu OPERA','Meas ALBA','Location','Southeast')
grid on
fprintf('   Sum oct = %8.2e     int %8.2e \n',ps(2),pint(2))

% plot deca
figure(5)
set(gcf,'color','w')
set(gca,'fontsize',16)
plot(traj0v10(2,:),pnv10(:,1),'-b','linewidth',2);hold on
plot(traj0(2,:),pn(:,1),'-r','linewidth',2);hold off
xlabel('S (mm)')
ylabel('Decapole')
legend('Simu OPERA','Meas ALBA','Location','Southeast')
grid on
fprintf('   Sum deca= %8.2e     int %8.2e \n',ps(1),pint(1))

% figure(100)
% set(gcf,'color','w')
% set(gca,'fontsize',16)
% subplot(1,1,4)
% plot(traj0v10(2,:),pnv10(:,4),'-b','linewidth',2);hold on
% plot(traj0(2,:),pn(:,4),'-r','linewidth',2);hold off
% xlabel('S (mm)')
% ylabel('Gradient')


% pnv10=pn;
% traj0v10=traj0;
% save dipole_v10.mat traj0v10 pnv10
return
 
figure(100); 
plot(xx,polyval(ps,xx*1e-3),'ob');hold on 
plot(xx,polyval(pint,xx*1e-3),'or');hold off 
grid on




