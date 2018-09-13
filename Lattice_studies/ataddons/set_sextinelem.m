function newring=set_sextinelem(ring,elemfam,scale)
% set sextupole component in elemfam 
% according to it qud strength K and disp value
% S=-K/D/2 in principle
newring=ring;
ind=findcells(newring,'FamName',elemfam);
% get dispersion
[TD] = twissring(newring, 0, ind,'chrom', 1e-8);
disp1  = cat(2,TD.Dispersion);
[TD] = twissring(newring, 0, ind+1,'chrom', 1e-8);
disp2  = cat(2,TD.Dispersion);
disp=(disp1(1,:)+disp2(1,:))/2;

for j=1:length(ind)
     K= newring{ind(j)}.PolynomB(2);
     newring{ind(j)}.PolynomB(3)= K*scale/disp(j)/2;
end