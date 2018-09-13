function [X] = atdata_full(RING)
% Some optics data for atmatch
% -I loc
% n1=13;n2=45; % 7BA
% n3=13+58;n4=43+58; % 6BA
% ns1=30;ns2=58+1;ns3=87;

% -I loc
n1=15;n2=48; % 7BA + drift
n3=15+62;n4=44+62; % 6BA +drift
% alpha=0 loc
ns1=32;ns2=62+1;ns3=93;

%
[TD, tunes, chrom] = twissring(RING, 0, 1:length(RING)+1,'chrom', 1e-8);
beta  = cat(1, TD.beta);
alpha = cat(1, TD.alpha);
%gamma = (1 + alpha.^2) ./ beta;
mu    = cat(1, TD.mu)/2/pi;
disp  = cat(2,TD.Dispersion);
[Emittance] = atemittance(RING,beta, alpha, disp);
%
% Data to fit
dmux1=abs(mu(n2,1)-mu(n1,1));
dmuz1=abs(mu(n2,2)-mu(n1,2));
dmux2=abs(mu(n3,1)-mu(n4,1));
dmuz2=abs(mu(n3,2)-mu(n4,2));

X=[tunes chrom beta(1,1)  beta(1,2)  disp(1,1) ...
   dmux1 dmuz1 dmux2 dmuz2 disp(1,15) alpha(ns1,1) alpha(ns1,2) ... 
  alpha(ns2,1) alpha(ns2,2) alpha(ns3,1) alpha(ns3,2) Emittance*1e12];

end
 
