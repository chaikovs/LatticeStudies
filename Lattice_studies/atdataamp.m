function [X] = atdataamp(RING)
% Some optics data for atmatch

[xdx, xdz, zdz]=nuamp_analyt_alex(RING);
X=[xdx xdz zdz]*1000;

end

