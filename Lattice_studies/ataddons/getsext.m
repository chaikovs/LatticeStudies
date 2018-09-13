function sext=getsext(ring)

ind=findcells(ring,'FamName','SXD1E');
sext(1)=ring{ind(1)}.PolynomB(3);
ind=findcells(ring,'FamName','SXF2E');
sext(2)=ring{ind(1)}.PolynomB(3);
ind=findcells(ring,'FamName','SXF1E');
sext(3)=ring{ind(1)}.PolynomB(3);
ind=findcells(ring,'FamName','SXF3E');
sext(4)=ring{ind(1)}.PolynomB(3);
ind=findcells(ring,'FamName','SXD2E');
sext(5)=ring{ind(1)}.PolynomB(3);



