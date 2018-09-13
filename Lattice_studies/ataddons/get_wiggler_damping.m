% based on Course A4: Damping Ring Design and Physics Issues Andy Wolsky
%
% some data
global GLOBVAL
E0   = GLOBVAL.E0;
Gamma=(E0+0.511e6)/0.511e6;
Brho = E0*1e-9/0.3;
Cq   = 3.8319e-13;
Cg   = 8.846e-5;

%
nw = 1.0;     % Wiggler number
Lw = 2.0;     % Wiggler length m
Bw = 0.92;    % Wiggler peak field T
rw = Brho/Bw;
lw = 0.020;  % Wiggler period m
kw = 2*pi/lw;
Dw = 0e-2;   % disersion at Wiggler loc m
bw = 2;      % min beta at wiggler centre m
bm = (bw*Lw/2 + (Lw/2)^3/3/bw)*2/Lw; % mean beta along Wiggler
%
I1w = -Lw/(2*(rw*kw)^2)*nw;
I2w =  Lw/2/rw^2*nw;
I3w =  4/3*Lw/rw^3*nw;
I4w =  0;
I5w =  4/(15*pi)*Lw*bm/kw^2/rw^5*nw + I3w*Dw^2/bw*nw;

%load('../Lattice-16fold/centre_opt1_7BA_type1_fit5', 'RING')
%load('../Lattice-16fold/centre_opt1_7BA_type1_fit9_atmatch3', 'RING')
%load Lattice-20fold/centre-20c-beam1-fit3-handfit1-beam4-beta4 
load('Lattice-save/SOLEIL_U_v9','RING')

par=atsummary(RING,'NoDisplay');
I1=par.integrals(1);
I2=par.integrals(2);
I3=par.integrals(3);
I4=par.integrals(4);
I5=par.integrals(5);

%
U0   = Cg*(E0*1e-9)^4*I2/2/pi;
eps0 = Cq*Gamma^2*I5/(I2-I4);
sige0= sqrt(Cq*Gamma^2*I3/(2*I2+I4));
Jx0  = 1 - I4/I2;
epsK0= eps0/(1+1/Jx0);
%
Uw   = Cg*(E0*1e-9)^4*I2w/2/pi;
epsw = Cq*Gamma^2*I5w/(I2w-I4w);
sigew= sqrt(Cq*Gamma^2*I3w/(2*I2w+I4w));
Jxw  = 1 - I4w/I2w;
epsKw= epsw/(1+1/Jxw);
%
U    = Cg*(E0*1e-9)^4*(I2+I2w)/2/pi;
eps  = Cq*Gamma^2*(I5+I5w)/(I2+I2w-I4-I4w);
sige = sqrt(Cq*Gamma^2*(I3+I3w)/(2*(I2+I2w)+I4+I4w));
Jx   = 1 - (I4+I4w)/(I2+I2w);
epsK = eps/(1+1/Jx);
%
fprintf('**********\n');
fprintf('                      U0         Eps       dE/E      Jx       EpsK \n');
fprintf('Bare lattice       : %8.2d   %8.2d  %8.2d  %8.2d  %8.2d\n',U0, eps0,sige0,Jx0,epsK0);
fprintf('Wiggler only       : %8.2d   %8.2d  %8.2d  %8.2d  %8.2d\n',Uw, epsw,sigew,Jxw,epsKw);
fprintf('lattice + wiggler  : %8.2d   %8.2d  %8.2d  %8.2d  %8.2d\n',U,  eps,sige,Jx,epsK);










