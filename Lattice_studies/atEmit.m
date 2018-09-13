function [X] = atEmit(RING,thetatot)
% Emittance and bending angle sum

% Set B63HE to get thetatot at
ind=findcells(RING,'FamName','B61HE');theta=RING{ind(1)}.BendingAngle*length(ind);
ind=findcells(RING,'FamName','B62HE');theta=theta+RING{ind(1)}.BendingAngle*length(ind);
ind=findcells(RING,'FamName','BQ62E');theta=theta+RING{ind(1)}.BendingAngle*length(ind);
ind=findcells(RING,'FamName','BQ63E');theta=theta+RING{ind(1)}.BendingAngle*length(ind);
ind=findcells(RING,'FamName','QF64E');theta=theta+RING{ind(1)}.BendingAngle*length(ind);
ind=findcells(RING,'FamName','B63HE');t3=(thetatot-theta)/length(ind);
for k=ind;RING{k}.BendingAngle=t3;end  
%
[TD] = twissring(RING, 0, 1:length(RING)+1,'chrom',1e-8);
beta  = cat(1, TD.beta);
alpha = cat(1, TD.alpha);
disp  = cat(2,TD.Dispersion);
[Emittance] = atemittance(RING,beta, alpha, disp);
%
X=[Emittance*1e12 disp(1,1)];

end
 
