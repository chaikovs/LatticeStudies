function ring=set_skewQinsext(ring,sext,ratio)

ind=findcells(ring,'FamName',sext);
for j=1:length(ind)
    strength=ring{ind(j)}.PolynomB(3);
    ring{ind(j)}.PolynomA(2)=strength*ratio;
end





