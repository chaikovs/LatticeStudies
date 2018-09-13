function track
%TRACK Summary of this function goes here
%   track orbit X over one thomX dipole
%   assume infinit radial size :  Bz only on axe
global E S Bz Bsp

E=2750;          % E MeV
brho=9.1747;
c=2.99792458e8;
beta=sqrt(1-(0.511/(E+0.511))^2);
%
Lbore=1.054;         % in mm
teta = pi/16;
% radia model
[S ,Bz, Bsp] = getEXCELfield1D('field/Bprofil-ans.dat',brho,teta,Lbore);
%[S ,Bz, Bsp] = getRADIAfield1D('field/Bprofil-s-ThomX-Dpole1.dat',brho,teta,Lbore);
%[S ,Bz, Bsp] = getMODELfield1D('field/Bprofil-s-ThomX-Dpole1.dat',brho,teta,Lbore);
len=length(S);

% Define initial conditions.
s1 = S(138);
s2 = S(275);
step=1000;
s=(s1:(s2-s1)/step:s2);
%
t=s/c;
X0 = [0 0.0001 s1  0  0.*c  c*beta]';

% Simulate the differential equation.
[t,X] = ode23('newton',t,X0);
X(step,4)

% basic matices track vertical
[Z] = matrix_coin(s,X0);

delta=X(:,4)./X(:,6); % xp
figure(10)
set(gca,'Fontsize',14)
plot(X(:,3)*1000,delta*1000)
xlabel('s (mm)')
ylabel('xp (mrad)')
xlim([s1 s2]*1000)
grid on



%
figure(2)
set(gca,'Fontsize',14)
plot(X(:,3)*1000,X(:,1)*1000)
xlim([s1 s2]*1000)
xlabel('s (mm)')
ylabel('x (mm)')
grid on

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
