function ring=setsext(ring,sext)

ind=findcells(ring,'FamName','SXD1E');
for j=1:length(ind)
    ring{ind(j)}.PolynomB(3)= sext(1);
end
ind=findcells(ring,'FamName','SXF2E');
for j=1:length(ind)
    ring{ind(j)}.PolynomB(3)= sext(2);
end
ind=findcells(ring,'FamName','SXF1E');
for j=1:length(ind)
    ring{ind(j)}.PolynomB(3)= sext(3);
end
ind=findcells(ring,'FamName','SXF3E');
for j=1:length(ind)
    ring{ind(j)}.PolynomB(3)= sext(4);
end
ind=findcells(ring,'FamName','SXD2E');
for j=1:length(ind)
    ring{ind(j)}.PolynomB(3)= sext(5);
end



