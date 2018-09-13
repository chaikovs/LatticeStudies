%function track_OPERA
%TRACK Summary of this function goes here
%   track orbit X over one thomX dipole
%   assume infinit radial size :  Bz only on axe
global E Xb Sb Bz 

E=70;          % E MeV
c=2.99792458e8;
beta=sqrt(1-(0.511/(E+0.511))^2);
%
radius=0.352; %m
Lbore=0.2765;         % in mm
teta = pi/8;
% radia model
[Xb ,Sb , Bz] = getOPERAfield2D('field/champ_dipole_OPERA_v10.table');
Bz=Bz*1.022;

[ traj0 , ds ] = make_traj(radius,0,teta);
ll0=length(traj0);
st0=traj0(2,ll0);

% Define initial conditions and ajust step to traj0
s1 = 0;
x1 = -0.00233
step=ll0-1;
s=(0:step)*ds;
%
t=s/c;
X0 = [x1 0.000 s1  0  0  c*beta]';

% Simulate the differential equation.
[t,X] = ode23('newton',t,X0);
tetaf=atan(X(:,4)./X(:,6))*180/pi; % teta
ll=length(tetaf);


tetaf0=90+atan((traj0(2,ll0)-traj0(2,ll0-1))./(traj0(1,ll0)-traj0(1,ll0-1)))*180/pi; % teta0


fprintf('Deviation = %f  / %f \n',tetaf0,-tetaf(ll))
%
figure(2)
set(gca,'Fontsize',14)
plot(traj0(2,:),traj0(1,:),'-k'); hold on
plot(X(:,3)*1000,X(:,1)*1000,'-b'); hold off
xlim([s1 st0])
xlabel('s (mm)')
ylabel('x (mm)')
legend('Ideal path','Tracked path')
grid on

%
dx=(X(:,1)*1000-traj0(1,:)');
figure(3)
set(gca,'Fontsize',14)
plot(X(:,3)*1000,dx)
xlim([s1 st0])
xlabel('s (mm)')
ylabel('dx (mm)')
grid on
return



figure(1)
plot(X(:,3)*1000,X(:,1)*1000,'-K'); hold on
image(Sb(:,1),Xb(1,:),-Bz'*100); 
plot(X(:,3)*1000,X(:,1)*1000,'-K'); hold off
xlim([0 300]);
ylim([-100 30 ])
xlabel('s (mm)')
ylabel('x (mm)')
grid on
return

%
figure(3)
set(gca,'Fontsize',14)
plot(X(:,3)*1000,X(:,2)*1000);hold on
plot(s*1000,Z(1,:)*1000,'-r');hold off
xlim([s1 s2]*1000)
xlabel('s (mm)')
ylabel('z (mm)');
legend('Tracking','Kick approximation K_1=0.17')
grid on
