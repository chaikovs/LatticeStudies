function [sposbend, Emittance_perbend]=atemittance_perbend(THERING,beta, alpha, disp)
%  emittance contribution per bend
global GLOBVAL
sum.e0 = GLOBVAL.E0*1e-9;
sum.gamma = sum.e0 / 0.51099906e-3;

% Synchrotron integral calculation
sum.integrals = zeros(1,6);

ii =0;
for i = 1:length(THERING),
    if isfield(THERING{i}, 'BendingAngle') && isfield(THERING{i}, 'EntranceAngle')
        ii = ii +1;
        rho = THERING{i}.Length/THERING{i}.BendingAngle;
        [dI1,dI2,dI3,dI4,dI5,curHavg1(ii), Dxavg(ii)] = ...
            calcRadInt(rho,THERING{i}.BendingAngle, ...
            alpha(i,1),beta(i,1),disp(1,i),disp(2,i),...
            THERING{i}.K,THERING{i}.EntranceAngle,THERING{i}.ExitAngle);
        
        sum.integrals(1) = sum.integrals(1) + dI1;
        sum.integrals(2) = sum.integrals(2) + dI2;
        sum.integrals(3) = sum.integrals(3) + dI3;
        % For general wedge magnets
        sum.integrals(4) = sum.integrals(4) + dI4;
        sum.integrals(5) = sum.integrals(5) + dI5;
        perbend(ii)=dI5;
        sposbend(ii)=findspos(THERING,i);
    end
end

% 
damping = 1 - sum.integrals(4)/sum.integrals(2);
% Emittance = 3.8319e-13*(sum.e0*1e3/0.510999).^2*sum.integrals(5)/(damping(1)*sum.integrals(2));
Emittance_perbend=3.8319e-13*(sum.e0*1e3/0.510999).^2*perbend./(damping(1)*sum.integrals(2));
%EnergySpread = sqrt(3.8319e-13*sum.gamma.^2*sum.integrals(3)/(2*sum.integrals(2) + sum.integrals(4)));



function [dI1,dI2,dI3,dI4,dI5,curHavg, Dxavg] = calcRadInt(rho,theta, a0,b0,D0,D0p,K1,th1,th2)
%[dI1,dI2,dI3,dI4,dI5,curHavg] = calcRadInt(rho,theta, a0,b0,D0,D0p,K1)
%calculate the contribution to the radiation integrals of a dipole.
%  INPUTS
%  rho, theta, radius and angle of the dipole
%  a0, b0, are horizontal alpha and beta at the entrance of the dipole,
%  D0, D0p are dispersion at the entrance of the dipole
%  K1, the gradient parameter in AT's convention, i.e., positive for
%  horizontal focusing, K1=0 by default
%  th1, th2, the entrance and exit angle, respectively, th1=th2= 0 [theta/2] by
%  default.
%

% If not combine dipole
if nargin == 6
    K1=0;
end

%
if nargin<9
    th1 = 0; %theta/2.0;
    th2 = 0; %theta/2.0;
end

% Edge focusing
M21 = tan(th1)/rho;
D0p = M21*D0+D0p;
a0 = -M21*b0+a0;

% split the dipole in N pieces
N = 100;
th = (0:N)/N*theta;

% Compute Twiss parameters inside dipole
for ii=1:length(th),
    [Dx(ii), Dxp(ii)] = calcdisp(rho, th(ii), D0, D0p, K1);
    [ax, bx] = calctwiss(rho, th(ii), a0, b0, K1);
    curHavg1(ii) = (Dx(ii)^2+(ax*Dx(ii)+bx*Dxp(ii))^2)/bx;
end

% Edge focusing
M21 = tan(th2)/rho;
Dxp(end) =  M21*Dx(end)+Dxp(end);
ax  = -M21*bx+ax;
curHavg1(end) = (Dx(end)^2+(ax*Dx(end)+bx*Dxp(end))^2)/bx;

% Average data
curHavg = ((curHavg1(1)+curHavg1(end))/2.0 + sum(curHavg1(2:end-1)))/(length(th)-1);
Dxavg   = ((Dx(1)+Dx(end))/2.0 + sum(Dx(2:end-1)))/(length(th)-1);

dI1 = ((Dx(1) + Dx(end))/2.0 + sum(Dx(2:end-1)))*theta/N;
dI2 = abs(theta/rho);
dI3 = abs(theta/rho^2);
dI4 = (1/rho^2 + 2*K1)*dI1  - (Dx(1)/rho^2*tan(th1) + Dx(end)/rho^2*tan(th2));
dI5 = curHavg*abs(theta/rho^2);

function [Dx, Dxp] = calcdisp(rho, theta, D0, D0p, K1)
%calcdisp - calculate dispersion function inside a combined-function dipole
%  INPUTS
%  1. rho - curvature radius
%  2. theta - angle
%  3. D0 - Horizontal dispersion function at the entrance
%  4. D0p - DErivative of  Horizontal dispersion function at the entrance
%  5. K1 - Focusing
%
% Transfert matrix of A wedge dipole p58 Handbook af Accelerator Physics
s = rho*theta;
if K1>-1/rho^2; %horizontal focusing
    sqK = sqrt(1/rho^2+K1);
    Dx  =  D0*cos(sqK*s) + D0p/sqK*sin(sqK*s)+(1-cos(sqK*s))/rho/sqK^2;
    Dxp = -D0*sqK*sin(sqK*s)+D0p*cos(sqK*s)+sin(sqK*s)/rho/sqK;
else %horizontal defocusing
    sqK = sqrt(-(1/rho^2+K1));
    Dx =  D0*cosh(sqK*s) + D0p/sqK*sinh(sqK*s)+(-1+cosh(sqK*s))/rho/sqK^2;
    Dxp = D0*sqK*sinh(sqK*s)+D0p*cosh(sqK*s)+sinh(sqK*s)/rho/sqK;
    
end

function [ax, bx] = calctwiss(rho, theta, a0, b0, K1)
% calctwiss calculate twiss function inside a combined-function dipole manget
%  INPUTS
%  1. rho - curvature radius
%  2. theta - angle
%  3. a0 - Horizontal alpha function at the entrance
%  4. b0 - Horizontal beta function at the entrance
%  5. K1 - Focusing
%
%  [beta ] = [  Mx11^2        -2*MX11*Mx12         Mx12^2   ] [beta0 ]
%  [alpha] = [ -Mx11*Mx21 Mx11*Mx22 + Mx11*Mx21   -Mx12*Mx22] [alpha0]
%  [gamma] = [  Mx21^2        -2*MX21*Mx22         Mx22^2   ] [gamma0]

Mx = calcMx(rho, K1,theta);
g0 = (1+a0^2)/b0;
twx2 = [Mx(1,1)^2, -2*Mx(1,1)*Mx(1,2), Mx(1,2)^2;
    -Mx(1,1)*Mx(2,1), Mx(1,1)*Mx(2,2)+Mx(1,2)*Mx(2,1),-Mx(1,2)*Mx(2,2);
    Mx(2,1)^2, -2*Mx(2,1)*Mx(2,2),Mx(2,2)^2] * [b0, a0, g0]';
ax = twx2(2);
bx = twx2(1);

function Mx = calcMx(rho,K1,theta)
% calcMx calculate transfer matrice of a combined-function dipole manget

s = rho*theta;
if K1>-1/rho^2; %horizontal focusing
    sqK = sqrt(1/rho^2+K1);
    Mx = [cos(sqK*s), sin(sqK*s)/sqK; -sqK*sin(sqK*s), cos(sqK*s)];
else %horizontal defocusing
    sqK = sqrt(-(1/rho^2+K1));
    Mx = [cosh(sqK*s), sinh(sqK*s)/sqK; sqK*sinh(sqK*s), cosh(sqK*s)];
end


