function p = seedquads(quaderr)

ao=getao;

quadfam={'QP1' 'QP2' 'QP3' 'QP4' 'QP31' 'QP41'};

for ii=1:length(quadfam)
   for jj=1:length(getfamilydata(quadfam{ii},'ElementList'));
    k0=getsp(quadfam{ii},elem2dev(quadfam{ii},jj),'physics');
    setsp(quadfam{ii},(k0+k0*quaderr*randn),elem2dev(quadfam{ii},jj),'physics');
    k=getsp(quadfam{ii},elem2dev(quadfam{ii},jj),'physics');
    %disp([quadfam{ii} ' element ' num2str(jj) ' value ' num2str(k0)])
   end
end
end

%plotbeta