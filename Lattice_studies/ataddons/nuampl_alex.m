function [nux,nuz]=nuampl_alex(ring,ampl,xz,dp,pl)
%ATNUAMPL	computes tune shift with amplitude
%[NUX,NUZ]=ATNUAMPL(RING,AMPLITUDE)
%
%	Computes tunes for the specified horizontal amplitudes
%
%[NUX,NUZ]=ATNUAMPL(RING,AMPLITUDE,3)
%
%	Computes tunes for the specified vertical amplitudes
% e.g.  atnuampl(esrf,0:.0002:0.01)
damp=1e-5;
if nargin < 3, xz=1; end
siza=size(ampl);
nampl=prod(siza);
p0=repmat([damp;0;damp;0;dp;0], 1,nampl);
p0(xz,:)=p0(xz,:)+ampl(:)';
p1=ringpass(ring,p0,128*4);

nux=reshape(findtune_multi(reshape(p1(1,:),nampl,[])',3),siza);
nuz=reshape(findtune_multi(reshape(p1(3,:),nampl,[])',3),siza);

if(pl)
    plot((ampl.*ampl)',[nux-nux(1);nuz-nuz(1)]','o-');
    legend('\nu_x','\nu_z');
    grid on
end