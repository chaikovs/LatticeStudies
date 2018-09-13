function [spos, hjklm]=calc_hjklm_alex(ring,j,k,l,m)
%dnudx=nuamp_analyt(ring,sym)
% add sextupole length
% return a vectot along the lattice

if nargin==1
   sym=1;
end

[lindat]=atlinopt(ring,0,1:length(ring)+1);

betaxy = cat(1, lindat.beta);
betax = betaxy(:,1);
betay = betaxy(:,2);
mu = cat(1, lindat.mu);
mux = mu(:,1);
muy = mu(:,2);
spos=cat(1,lindat.SPos);
%now we need to add up h_jklm contribution from the sextupoles.
% 
%basically, we define
%hjklm = sum_u b3L_u beta_x^(j+k)/2 beta_y^(l+m)/2
%e^(i(j-k)phi_x+i(l-m)phi_y)
 
index=findcells(ring,'PolynomB');
hjklm=ones(length(index)+1,1);
hjklm(1)=0;
for q=1:length(index)
    b3l=ring{index(q)}.PolynomB(3)*ring{index(q)}.Length;
    hjklm(q+1)=hjklm(q)+b3l*betax(index(q))^((j+k)/2)*betay(index(q))^((j+k)/2)*...
        exp(1i*((j-k)*mux(index(q))+(l-m)*muy(index(q))));
end
hjklm(q+2)=hjklm(q+1);
spos=[0 ; spos(index) ; spos(length(spos))];





