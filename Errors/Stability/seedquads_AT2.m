function rerr = seedquads_AT2(ring, quaderrsigma)

% ao=getao;
% 
% quadfam={'QP1' 'QP2' 'QP3' 'QP4' 'QP31' 'QP41'};
% 
% for ii=1:length(quadfam)
%    for jj=1:length(getfamilydata(quadfam{ii},'ElementList'));
%     k0=getsp(quadfam{ii},elem2dev(quadfam{ii},jj),'physics');
%     setsp(quadfam{ii},(k0+k0*quaderr*randn),elem2dev(quadfam{ii},jj),'physics');
%     k=getsp(quadfam{ii},elem2dev(quadfam{ii},jj),'physics');
%     %disp([quadfam{ii} ' element ' num2str(jj) ' value ' num2str(k0)])
%    end
% end
% end

%plotbeta



% get indexes
indm=find(atgetcells(ring,'FamName','BPMx'));
indq=find(atgetcells(ring,'Class','Quadrupole'));
indd=find(atgetcells(ring,'Class','Bend'));

%BareRING = atsetfieldvalues(ring,findcells(ring,'FamName','QP1'), 'PolynomB',{1,2},0);

rerr=atsetrandomerrors(...
    ring,...
    indq,...
    indm,...
    1,...
    quaderrsigma,...
    2,...
    'dpb2');

end