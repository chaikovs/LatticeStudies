function [RING] = AB_bend(RING,r1,r2,type)
%UNTITLED Summary of this function goes here
% symmetrize the to AB quad component

if strcmp(type,'7BA')

elseif strcmp(type,'6BA')
    ind1=findcells(RING,'FamName','BQ62');
    T1=RING{ind1(1)}.BendingAngle;
    ind2=findcells(RING,'FamName','BQ63');
    T2=RING{ind2(1)}.BendingAngle;
end
%
TM=(T1+T2)/2;
%
T1=T1+r1*(TM-T1);
T2=T2+r2*(TM-T2);
%

for i=1:length(ind1);RING{ind1(i)}.BendingAngle=T1;end
for i=1:length(ind2);RING{ind2(i)}.BendingAngle=T2;end

end

