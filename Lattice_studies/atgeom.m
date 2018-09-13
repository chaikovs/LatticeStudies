function [X] = atgeom(RING)
% Some optics data for atmatch
% -I loc
ind=findcells(RING,'BendingAngle');
theta=0;
for i=ind
    theta=theta+RING{i}.BendingAngle;
end

X=theta;
end

