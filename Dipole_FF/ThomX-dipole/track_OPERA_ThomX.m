%function track_OPERA_ThomX
%TRACK Summary of this function goes here
%   track orbit X over one thomX dipole
%   assume infinit radial size :  Bz only on axe
global E Xb Sb Bz 

E=50;          % E MeV
c=2.99792458e8;
beta=sqrt(1-(0.511/(E+0.511))^2);
%
radius=0.352; %m
Lbore=0.2765;         % in mm
teta = pi/8;
% radia model
%[Xb ,Sb , Bz] = getOPERAfield2D('field/champ_dipole_OPERA_v31.table');
% [Xb ,Sb , Bz] = getOPERAfield2D('field/champ_dipole_OPERA_v33bis.table');
% [Xb ,Sb , Bz] = getOPERAfield2D('field/dipole33_XC10.table');
% Bz=Bz*0.7965;

[Xb,Sb,Bz] = getMEASUREDfield2D('field/20160208a_THOMX#009_fieldmap_160A.xlsx');
Bz=Bz*0.996; % *0.996 pour ~50 MeV 160A

% tetapole=pi/9;  
% [pole] = make_pole(radius,0.050,2*tetapole);
% figure(11)
% set(gcf,'color','w')
% set(gca,'fontsize',18)
% plot(pole(2,:),pole(1,:),'-w','Linewidth',3);
% hold on
% imagesc(Sb(:,1),Xb(1,:),Bz');
% shading interp;
% xlim([-250 250])
% ylim([-100 30])
% xlabel('s (mm)')
% ylabel('x (mm)')
% grid on
% hold off
% colorbar

[ traj0 , ds ] = make_traj(radius,0,teta);
ll0=length(traj0);
st0=traj0(2,ll0);

% Define initial conditions and ajust step to traj0
s1 = 0;
% x1 = -0.0027 + 0.0000296;
x1 = 0;%-0.0027 + 0.0000296;
step=ll0-1;
s=(0:step).*ds(1);

%
t=s/c;
X0 = [x1 0.000 s1  0  0  c*beta]';

% Simulate the differential equation.
[t,X] = ode23('newton',t,X0);

BB= [];
for i=1:length(X)
    BB = [BB; getfield2(X(i,:))];
end;    
sss = sqrt(diff(X(:,1)).^2 + diff(X(:,3)).^2);
sssn = [0;sss];


tetaf=atan(X(:,4)./X(:,6))*180/pi; % teta
ll=length(tetaf);


tetaf0=90+atan((traj0(2,ll0)-traj0(2,ll0-1))./(traj0(1,ll0)-traj0(1,ll0-1)))*180/pi; % teta0
%90+atan((S0(end)-S0(end-1)).*1e3./(traj0(1,ll0)-traj0(1,ll0-1)))*180/pi

fprintf('Deviation = %f  / %f \n',tetaf0,-tetaf(ll))
%
figure(2)
set(gcf,'color','w')
set(gca,'fontsize',16)
plot(traj0(2,:),traj0(1,:),'-k'); hold on
plot(X(:,3)*1000,X(:,1)*1000,'-b'); hold off
%xlim([s1 st0])
xlabel('s (mm)')
ylabel('x (mm)')
legend('Ideal path','Tracked path')
grid on
print('dipole_track_ideal_path','-dpng','-r300')
%
dx=(X(:,1)*1000-traj0(1,:)');
figure(3)
set(gcf,'color','w')
set(gca,'fontsize',16)
plot(X(:,3)*1000,dx)
xlim([s1 st0])
xlabel('s (mm)')
ylabel('Tracked - Ideal path (mm)')
grid on
print('dipole_track_ideal_path_DIFF','-dpng','-r300')

bz   =interp2(Xb,Sb,Bz,traj0(1,:)  ,traj0(2,:));

Yr=X(:,3)*1000;
X0=traj0(1,:)';
Y0=traj0(2,:)';
S0=cumsum(ds) - ds(1);
Xr=X(:,1)*1000;
Sr= cumsum(sssn)*1000;
dX=dx;
Bztrack = BB(:,2);
Bz0 = bz;

intbz0=trapz(Sr/1000,Bztrack) ;
gap   =0.042 ;

fprintf('  Peak field  %g  T \n',Bztrack(1) )
fprintf('  Integrated field  %g  mTm \n',intbz0 *2*1000)
fprintf('  Equivalent length %g  mm \n',intbz0/BB(1,2)*2*1000)
K = 1./(gap.*BB(1,2).^2) * trapz(Sr/1000, Bztrack.*(BB(1,2) - Bztrack));
fprintf('  Fringe coef K = %g  \n',K)

%save ThomX_dipole_traj160A.mat S Bztrack intbz0
%save ThomX_dipole_trajectoire.mat S X0 Xr dX
%save ThomX_dipole_trajectoryALL160A_m20mm.mat X0 Y0 S0 Xr Yr Sr Bztrack Bz0 dX



figure(10)
set(gcf,'color','w')
set(gca,'fontsize',16)
plot(cumsum(sssn)*1000, BB(:,2),'r.-')
hold on
plot((cumsum(ds)-ds(1))*1000,bz,'.-m')
%plot(cumsum(ds)*1000,bz,'.-k')
xlabel('s-traj (mm)')
ylabel('B (T)')
legend('Tracking','Ideal trajectory')
grid on
hold off
print('dipole_track_ideal_BS','-dpng','-r300')

figure(1)
set(gcf,'color','w')
set(gca,'fontsize',16)
plot(X(:,3)*1000,X(:,1)*1000,'-K'); hold on
image(Sb(:,1),Xb(1,:),Bz'*100); 
plot(X(:,3)*1000,X(:,1)*1000,'-K'); hold off
xlim([0 300]);
ylim([-100 30 ])
xlabel('s (mm)')
ylabel('x (mm)')
grid on
print('dipole_track_image','-dpng','-r300')

return

%
figure(4)
set(gca,'Fontsize',14)
plot(X(:,3)*1000,X(:,2)*1000);hold on
plot(s*1000,Z(1,:)*1000,'-r');hold off
xlim([s1 s2]*1000)
xlabel('s (mm)')
ylabel('z (mm)');
legend('Tracking','Kick approximation K_1=0.17')
grid on
