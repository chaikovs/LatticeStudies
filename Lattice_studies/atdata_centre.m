function [X] = atdata_centre(RING)
% Some optics data for atmatch

%
[TD, tunes, chrom] = twissring(RING, 0, 1:length(RING)+1,'chrom', 1e-8);
beta  = cat(1, TD.beta);
alpha = cat(1, TD.alpha);
%gamma = (1 + alpha.^2) ./ beta;
mu    = cat(1, TD.mu)/2/pi;
disp  = cat(2,TD.Dispersion);
[Emittance] = atemittance(RING,beta, alpha, disp);

% Data to fit
X=[tunes  beta(1,1)  beta(1,2) disp(1,1) Emittance*1e12];
%X={tunes(1), tunes(1), beta(1,1),  beta(1,2), disp(1,1), Emittance*1e12};

end

