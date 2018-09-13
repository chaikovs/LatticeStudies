function [newring]=set_sextlength(ring,lsext)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

newring=ring;
ind=findcells(newring,'Class','Sextupole');

for j=1:1:length(ind)
   lsext0=newring{ind(j)}.Length;
   newring{ind(j)}.Length=lsext;
   newring{ind(j)}.PolynomB(3)=newring{ind(j)}.PolynomB(3)*lsext0/lsext;
   k=0; while ~strcmp(newring{ind(j)+k}.Class,'Drift');k=k-1;end
   newring{ind(j)+k}.Length=newring{ind(j)+k}.Length-(lsext-lsext0)/2;
   k=0; while ~strcmp(newring{ind(j)+k}.Class,'Drift');k=k+1;end
   newring{ind(j)+k}.Length=newring{ind(j)+k}.Length-(lsext-lsext0)/2;
end
end

