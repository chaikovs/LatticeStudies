function [RING] = AB_neq(RING,r1,r2,type)
%UNTITLED Summary of this function goes here
% symmetrize the to AB quad component

if strcmp(type,'7BA')
    ind1=findcells(RING,'FamName','BQ62E');
    K1=RING{ind1(1)}.PolynomB(2);
    ind2=findcells(RING,'FamName','BQ63E');
    K2=RING{ind2(1)}.PolynomB(2);
elseif strcmp(type,'6BA')
    ind1=findcells(RING,'FamName','BQ62');
    K1=RING{ind1(1)}.PolynomB(2);
    ind2=findcells(RING,'FamName','BQ63');
    K2=RING{ind2(1)}.PolynomB(2);
end
%
%
K1=K1*r1;
K2=K2*r2;
%

for i=1:length(ind1);RING{ind1(i)}.PolynomB(2)=K1;end
for i=1:length(ind2);RING{ind2(i)}.PolynomB(2)=K2;end

end

