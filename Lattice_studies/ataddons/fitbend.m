function [ RING ] = fitbend(RING, thetatot, r, type )
%fit bend deviation to get Dx and theta total
%   for 6BA or 7BA 

if strcmp(type,'7BA')
    r1=r(1);
    r2=r(2);
    r3=r(3);
    r4=r(4); %antibend
    ind=findcells(RING,'FamName','B61HE');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r1;end
    theta=RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B62HE');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r2;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B63HE');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r3;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','BQ62E');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r4;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','BQ63E');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r4;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B64HE');
    t4=(thetatot-theta)/length(ind);
    for k=ind;RING{k}.BendingAngle=t4;end
    %theta=theta+RING{k}.BendingAngle*length(ind);
    
elseif strcmp(type,'7BAsym') % B64HE = B63HE
    
    r1=r(1);
    r2=r(2);
    r3=r(3);
    ind=findcells(RING,'FamName','B61HE');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r1;end
    theta=RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B62HE');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r2;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','BQ62E');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r3;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','BQ63E');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r3;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B63HE');
    t3=(thetatot-theta)/length(ind);
    for k=ind;RING{k}.BendingAngle=t3;end
    
  elseif strcmp(type,'7BAsym16') % B64HE = B63HE
    
    r1=r(1);
    r2=r(2);
    r3=r(3);
    r4=r(4);
    ind=findcells(RING,'FamName','B61HE');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r1;end
    theta=RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B62HE');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r2;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','BQ62E');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r3;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','BQ63E');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r3;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','QF64E');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r4;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B63HE');
    t3=(thetatot-theta)/length(ind);
    for k=ind;RING{k}.BendingAngle=t3;end  
    
  elseif strcmp(type,'7BAsym16QF65') % B64HE = B63HE + QF65E
    
    r1=r(1);
    r2=r(2);
    r3=r(3);
    r4=r(4);
    r5=r(5);
    ind=findcells(RING,'FamName','B61HE');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r1;end
    theta=RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B62HE');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r2;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','BQ62E');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r3;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','BQ63E');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r3;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','QF65E');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r4;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','QF64E');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r5;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B63HE');
    t3=(thetatot-theta)/length(ind);
    for k=ind;RING{k}.BendingAngle=t3;end  
    
elseif strcmp(type,'7BAsdl')
    r1=r(1);  % sdm side
    r2=r(2);  % sdl side
    r3=r(3);
    r4=r(4);
    r5=r(5);  %antibend
    ind=findcells(RING,'FamName','B61HE');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r1;end
    theta=RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B61HL');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r2;end  
    theta= theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B62HE');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r3;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B63HE');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r4;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','BQ62E');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r5;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','BQ63E');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r5;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B64HE');
    t4=(thetatot-theta)/length(ind);
    for k=ind;RING{k}.BendingAngle=t4;end
    %theta=theta+RING{k}.BendingAngle*length(ind);
    
 elseif strcmp(type,'7BAsymsdl') % B64HE = B63HE
    r1=r(1);% sdm side
    r2=r(2);% sdl side
    r3=r(3);
    r4=r(4);%antibend   
    ind=findcells(RING,'FamName','B61HE');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r1;end
    theta=RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B61HL');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r2;end  
    theta= theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B62HE');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r3;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','BQ62E');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r4;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','BQ63E');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r4;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B63HE');
    t3=(thetatot-theta)/length(ind);
    for k=ind;RING{k}.BendingAngle=t3;end
    
elseif strcmp(type,'6BA')
    
    r1=r(1);
    r2=r(2);
    r3=r(3);
    ind=findcells(RING,'FamName','B61H');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r1;end
    theta=RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B62H');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r2;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','BQ62');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r3;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','BQ63');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r3;end
    theta=theta+RING{k}.BendingAngle*length(ind);
    ind=findcells(RING,'FamName','B63H');
    t3=(thetatot-theta)/length(ind);
    for k=ind;RING{k}.BendingAngle=t3;end
    %theta=theta+RING{k}.BendingAngle*length(ind);
    
elseif strcmp(type,'8BA')   
    
    r1=r(1);
    r2=r(2);
    r3=r(3);
    ind=findcells(RING,'FamName','B64H');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r1;end
    theta=RING{k}.BendingAngle*length(ind)
    ind=findcells(RING,'FamName','B61H');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r2;end
    theta=theta+RING{k}.BendingAngle*length(ind)
    ind=findcells(RING,'FamName','B62H');
    for k=ind;RING{k}.BendingAngle=RING{k}.BendingAngle*r3;end
    theta=theta+RING{k}.BendingAngle*length(ind)
    t4=(thetatot-theta)/4;
    ind=findcells(RING,'FamName','B63H');
    for k=ind;RING{k}.BendingAngle=t4;end
    theta=theta+RING{k}.BendingAngle*length(ind) 
    
else
end

end

