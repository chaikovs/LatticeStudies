%Compute ring life-time
clear
fprintf(' ******** Lifetime ************ \n')
e0=0.511; % MeV
E=2750;  % MeV
gam=(E+e0)/e0;
betaxm=8;
betazm=8;
Q=1.e-9;     % Coulomb per bucket
Circ=354; % m;
epsx=20e-12; %Upgrade
epsz=20e-12;
% epsx=4e-9;  % Actual
% epsz=50e-12;
acc= 0.01; % energy acceptance

% Vacuum lifetime
r0=2.82e-15;
c=3e8;
avo=6.6e22;
Z=7;

a=0.010;
b=0.008;
p=1e-9; % mbar
rho=avo*p;
term=2*pi*r0^2*Z^2*rho*c/gam^2*(betaxm*betaxm/a^2 + betazm*betazm/b^2);
tau=1/term/3600;
fprintf('Lifetime gas = %5.1d h (%5.0f s) at %d MeV\n',tau,tau*3600,E)


% Touschek Lifetime 
e=1.6e-19;
Np=Q/e;
%sigl=2.3e-3 ;  % en m
sigl=0.8;  % en m
sigx=sqrt(epsx*betaxm);
sigz=sqrt(epsz*betazm);
sigxp=epsx/sigx;
x=(acc/sigxp/gam)^2;
F=intF(x);

tau=8*pi*sigl*sigx*sigz*sigxp*acc^2*gam^3/Np/r0^2/c/F/3600;
fprintf('Lifetime Tou = %5.1d h (%5.0f s) at %d MeV  and  %d  nC \n',tau,tau*3600,E,Np*e*1e9)





