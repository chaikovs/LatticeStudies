
%%
ring = thomx_lattice();

indDip=find(atgetcells(ring,'Class','Bend'));

PB_dip = atgetfieldvalues(ring,indDip,'PolynomB',{1,1})

% 0.296233/0.785398163397448
%%

r_optic=atsetfieldvalues(ring, indDip,'PolynomB',{1,1}, 0.37717557);

%%


