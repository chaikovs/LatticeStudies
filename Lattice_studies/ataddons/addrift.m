function [RING] = addrift(RING,ds,type)
% add octupoles at SX location
%   for 6BA or 7BA 

if strcmp(type,'7BA')
    
    D=atelem('drift','FamName','D621','Length',ds);
    sx=findcells(RING,'FamName','SXD1E');
    lr=length(RING);RING=[RING(1:sx(1)-1) D RING(sx(1) : sx(2)) D  RING(sx(2)+1:lr)];
    sd=findcells(RING,'FamName','D62');
    for i=1:length(sd);RING{sd(i)}.Length=RING{sd(i)}.Length-ds;end
  
    D=atelem('drift','FamName','D641','Length',ds);
    sx=findcells(RING,'FamName','SXD2E');
    lr=length(RING);RING=[RING(1:sx(1)) D RING(sx(1)+1 : sx(2)-1) D  RING(sx(2):lr)];
    sd=findcells(RING,'FamName','D64');
    for i=1:length(sd);RING{sd(i)}.Length=RING{sd(i)}.Length-ds;end
        
elseif strcmp(type,'6BA')
    
    D=atelem('drift','FamName','D621','Length',ds);
    sx=findcells(RING,'FamName','SXD1');
    lr=length(RING);RING=[RING(1:sx(1)-1); D; RING(sx(1) : sx(2)); D;  RING(sx(2)+1:lr)];
    sd=findcells(RING,'FamName','D62');
    for i=1:length(sd);RING{sd(i)}.Length=RING{sd(i)}.Length-ds;end
  
    D=atelem('drift','FamName','D641','Length',ds);
    sx=findcells(RING,'FamName','SXD2');
    lr=length(RING);RING=[RING(1:sx(1)); D; RING(sx(1)+1 : sx(2)-1); D;  RING(sx(2):lr)];
    sd=findcells(RING,'FamName','D64');
    for i=1:length(sd);RING{sd(i)}.Length=RING{sd(i)}.Length-ds;end  
    
else
end
end

