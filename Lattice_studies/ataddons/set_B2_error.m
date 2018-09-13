function [ring]=set_B2_error(ring,class,dksk)
%UNTITLED Summary of this function goes here
%   Apply a relative quad error to element class : Quadrupole, Bend
%   Simple Gaussian distribution

ind=findcells(ring,'Class',class);
for j=1:length(ind)
    qtest=ring{ind(j)}.PolynomB(2);
    if qtest~=0
        ring{ind(j)}.PolynomB(2)= qtest*(1+randn(1)*dksk);
    end
end
end

