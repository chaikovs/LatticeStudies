function [xmax] = atdynap_fast(ring,xlist,zfix,dpp,nt)
%Compute the dynamic aperture for z fixed on negative side
%
if nargin < 5, nt=300; end
if nargin < 4, dpp=0.0; end

% if isnumeric(dpp)
% clorb=[findorbit4(ring,dpp);dpp;0];
% else
%    clorb=findorbit6(ring);
% end

xmax = 0.0;
for i=1:length(xlist)
   rin=[xlist(i);0;zfix;0;dpp;0];
   [dummy,lost]=ringpass(ring,rin,nt,'KeepLattice'); %#ok<ASGLU>
   if lost, break; end
end
if i>1 ; xmax=xlist(i-1);end
end

