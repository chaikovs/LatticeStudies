function [X] = atdata_thomx(RING)
% Some optics data for atmatch
% -I loc
ind=findcells(RING,'FamName','SD0');n1=ind(1);n2=ind(2);
%
[TD, tunes, chrom] = twissring(RING, 0, 1:length(RING)+1,'chrom', 1e-8);
beta  = cat(1, TD.beta);
alpha = cat(1, TD.alpha);
%gamma = (1 + alpha.^2) ./ beta;
mu    = cat(1, TD.mu)/2/pi;
disp  = cat(2,TD.Dispersion);
[Emittance, EnergySpread, Jx] = atemittance(RING,beta, alpha, disp);
%
% Data to fit
dmux=mu(n2,1)-mu(n1,1);
dmuz=mu(n2,2)-mu(n1,2);
X=[tunes chrom beta(1,1)  beta(1,2)  disp(1,1) dmux dmuz disp(1,n1) Emittance*1e12];

%MI=TD(n2).M44*inv(TD(n1).M44)

end

