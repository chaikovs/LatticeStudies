function oct=getoct(ring)

ind=findcells(ring,'FamName','OCTF1');
oct(1)=ring{ind(1)}.PolynomB(4);
ind=findcells(ring,'FamName','OCTD1');
oct(2)=ring{ind(1)}.PolynomB(4);
ind=findcells(ring,'FamName','OCTD2');
oct(3)=ring{ind(1)}.PolynomB(4);




