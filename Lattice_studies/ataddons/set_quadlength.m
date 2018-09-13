function [newring]=set_quadlength(ring,name,len)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

newring=ring;
ind=findcells(newring,'FamName',name);

for j=1:1:length(ind)
   len0=newring{ind(j)}.Length;
   newring{ind(j)}.Length=len;
   newring{ind(j)}.PolynomB(2)=newring{ind(j)}.PolynomB(2)*len0/len;
   k=0; while ~strcmp(newring{ind(j)+k}.Class,'Drift');k=k-1;end
   newring{ind(j)+k}.Length=newring{ind(j)+k}.Length-(len-len0)/2;
   k=0; while ~strcmp(newring{ind(j)+k}.Class,'Drift');k=k+1;end
   newring{ind(j)+k}.Length=newring{ind(j)+k}.Length-(len-len0)/2;
end
end

