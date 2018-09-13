function dx=getreversebend(ring)
% dx for BQ62E BQ63E QF64E

ind=findcells(ring,'FamName','BQ62E');
L=ring{ind(1)}.Length;
T=ring{ind(1)}.BendingAngle;
K=ring{ind(1)}.PolynomB(2);
dx(1)=T/K/L;

ind=findcells(ring,'FamName','BQ63E');
L=ring{ind(1)}.Length;
T=ring{ind(1)}.BendingAngle;
K=ring{ind(1)}.PolynomB(2);
dx(2)=T/K/L;

ind=findcells(ring,'FamName','QF64E');
L=ring{ind(1)}.Length;
T=ring{ind(1)}.BendingAngle;
K=ring{ind(1)}.PolynomB(2);
dx(3)=T/K/L;