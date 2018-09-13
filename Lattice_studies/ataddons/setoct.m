function ring=setoct(ring,oct)

ind=findcells(ring,'FamName','OCTF1');
for j=1:length(ind)
    ring{ind(j)}.PolynomB(4)= oct(1);
end
ind=findcells(ring,'FamName','OCTD1');
for j=1:length(ind)
    ring{ind(j)}.PolynomB(4)= oct(2);
end
ind=findcells(ring,'FamName','OCTD2');
for j=1:length(ind)
    ring{ind(j)}.PolynomB(4)= oct(3);
end




