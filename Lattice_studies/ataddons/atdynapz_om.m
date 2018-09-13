function [zmaxlist,dplist] = atdynapz_om(ring,xfix,zlist,dplist,nt)
%Compute the dynamic aperture for z fixed on negative side
%along momentum
%
if nargin < 5, nt=300; end

zmaxlist(1:length(dplist))=0;
for i=1:length(dplist)
    [zmaxlist(i)] =atdynapz_fast(ring,xfix,zlist,dplist(i),nt);
end

end

