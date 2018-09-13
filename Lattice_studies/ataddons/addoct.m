function [RING] = addoct(RING,oct,type)
% add octupoles at SX location
%   for 6BA or 7BA 

if strcmp(type,'7BA')
    
    SO=oct;
    order=4;
    polb(1,order)=SO(1);
    OCTF1=atelem('marker','FamName','OCTF1','PolynomA',zeros(1,order),...
        'PolynomB',polb,'MaxOrder',order-1);
    OCTF1.PassMethod='ThinMPolePass';
    polb(1,order)=SO(2);
    OCTD1=atelem('marker','FamName','OCTD1','PolynomA',zeros(1,order),...
        'PolynomB',polb,'MaxOrder',order-1);
    OCTD1.PassMethod='ThinMPolePass';
    polb(1,order)=SO(3);
    OCTD2=atelem('marker','FamName','OCTD2','PolynomA',zeros(1,order),...
        'PolynomB',polb,'MaxOrder',order-1);
    OCTD2.PassMethod='ThinMPolePass';
    
    sx=findcells(RING,'FamName','SXF1E');
    lr=length(RING);RING=[RING(1:sx(1)); OCTF1; RING(sx(1)+1:sx(2)); OCTF1;  RING(sx(2)+1:lr)];
    sx=findcells(RING,'FamName','SXD1E');
    lr=length(RING);RING=[RING(1:sx(1)); OCTD1; RING(sx(1)+1:sx(2)); OCTD1;  RING(sx(2)+1:lr)];
    sx=findcells(RING,'FamName','SXD2E');
    lr=length(RING);RING=[RING(1:sx(1)); OCTD2; RING(sx(1)+1:sx(2)); OCTD2;  RING(sx(2)+1:lr)];
    
elseif strcmp(type,'6BA')
    
    SO=oct;
    order=4;
    polb(1,order)=SO(1);
    OCTF1=atelem('marker','FamName','OCTF1','PolynomA',zeros(1,order),...
                   'PolynomB',polb,'MaxOrder',order-1);
    OCTF1.PassMethod='ThinMPolePass';
    polb(1,order)=SO(2);
    OCTD1=atelem('marker','FamName','OCTD1','PolynomA',zeros(1,order),...
                   'PolynomB',polb,'MaxOrder',order-1);
    OCTD1.PassMethod='ThinMPolePass';
    polb(1,order)=SO(3);
    OCTD2=atelem('marker','FamName','OCTD2','PolynomA',zeros(1,order),...
                   'PolynomB',polb,'MaxOrder',order-1);
    OCTD2.PassMethod='ThinMPolePass';
    
    sx=findcells(RING,'FamName','SXF1');
    lr=length(RING);RING=[RING(1:sx(1)) OCTF1 RING(sx(1)+1:sx(2)) OCTF1  RING(sx(2)+1:lr)];
    sx=findcells(RING,'FamName','SXD1');
    lr=length(RING);RING=[RING(1:sx(1)) OCTD1 RING(sx(1)+1:sx(2)) OCTD1  RING(sx(2)+1:lr)];
    sx=findcells(RING,'FamName','SXD2');
    lr=length(RING);RING=[RING(1:sx(1)) OCTD2 RING(sx(1)+1:sx(2)) OCTD2  RING(sx(2)+1:lr)];
    
else
end
end

