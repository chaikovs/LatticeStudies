function [twissdatain] = atdatain(RING,index)
%UNTITLED Summary of this function goes here
%  get entrance twiss data from cyclic sol
%

[TD, tunes] = twissring(RING, 0, 1:length(RING)+1,'chrom', 1e-8);
beta  = cat(1, TD.beta);
alpha = cat(1, TD.alpha);
mu    = cat(1, TD.mu)/2/pi;
disp  = cat(2,TD.Dispersion);
%

ind=1;
if nargin==2;ind=index;end
twissdatain.ElemIndex=ind;
twissdatain.SPos=0;
twissdatain.ClosedOrbit=zeros(4,1);
twissdatain.M44=eye(4);
twissdatain.beta= [beta(ind,1)  beta(ind,2)];
twissdatain.alpha=[alpha(ind,1)  alpha(ind,2)];
twissdatain.mu= [0 0];
twissdatain.Dispersion= [disp(1,ind) disp(2,ind) disp(3,ind) disp(4,ind)]';

end

