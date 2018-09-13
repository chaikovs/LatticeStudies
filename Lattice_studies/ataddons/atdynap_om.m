function [xmaxlist,dplist] = atdynap_om(ring,xlist,zfix,dplist,nt)
%Compute the dynamic aperture for z fixed on negative side
%along momentum
%
if nargin < 5, nt=300; end

xmaxlist(1:length(dplist))=0;
for i=1:length(dplist)
    [xmaxlist(i)]=atdynap_fast(ring,xlist,zfix,dplist(i),nt);
end

end

