% step tune

[TD, tunes, chrom] = twissring(RING, 0, 1:length(RING)+1,'chrom', 1e-8);
beta  = cat(1, TD.beta);

index=findcells(RING,'PolynomB');
li=length(index);
p=0;m=0;
for q=1:li
   S=RING{index(q)}.PolynomB(2);
   L=RING{index(q)}.Length;
   if S>=0 ; 
       p=p+1;
       KP(p) =S*L;
       betax(p)=beta(index(q),1);
   else
       m=m+1;
       KM(m)=S*L;
       betaz(m)=beta(index(q),1);
   end;
   %K(q)=S*L; % KL
end


T=[sum(betax.*KP)  sum(betax.*KM) ; -sum(betaz.*KP)  -sum(betaz.*KM)]