function [RING] = lengthbend(RING, r, type )
%fit bend deviation to get Dx and theta total
%   for 6BA or 7BA 

if strcmp(type,'7BA')

elseif strcmp(type,'6BA')
   
    ind=findcells(RING,'FamName','B61H');
    for k=ind;
        L=RING{k}.Length;
        RING{k}.Length=L*r;
        RING{k}.PolynomB(2)=RING{k}.PolynomB(2)/r;
        RING{k}.K=RING{k}.PolynomB(2);
    end
    ind=findcells(RING,'FamName','D65');
    for k=ind;
        RING{k}.Length=RING{k}.Length+L*(1-r);
    end
    ind=findcells(RING,'FamName','D62');
    for k=ind;
        RING{k}.Length=RING{k}.Length+L*(1-r);
    end    
   
else
end

end

