function [X] = atdatachro(RING)
% Some optics data for atmatch
% -I loc
n1=14;n2=43; % 6BA 3Q
%n1=10;n2=39; % 6BA  1Q
%n1=14;n2=45; % 7BA
%n1=10;n2=41; % 7BA 1Q
%n1=13;n2=43; % 6BA
%n1=15;n2=49; % 7BA OCT
%
RING=fitchrom_alex(RING,[0 0],'SXD2' ,'SXF1' );% 6BA
%
[TD, tunes] = twissring(RING, 0, 1:length(RING)+1,'chrom', 1e-8);
beta  = cat(1, TD.beta);
alpha = cat(1, TD.alpha);
%gamma = (1 + alpha.^2) ./ beta;
mu    = cat(1, TD.mu)/2/pi;
disp  = cat(2,TD.Dispersion);
[Emittance] = atemittance(RING,beta, alpha, disp);
%
% Data to fit
dmux=mu(n2,1)-mu(n1,1);
dmuz=mu(n2,2)-mu(n1,2);
X=[tunes beta(1,1)  beta(1,2)  disp(1,1) dmux dmuz disp(1,12) Emittance*1e12];

end

