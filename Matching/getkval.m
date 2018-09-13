
function [qall] = getkval(ring)
qp1ind     = findcells(ring,'FamName','QP1');
qp2ind     = findcells(ring,'FamName','QP2');
qp3ind     = findcells(ring,'FamName','QP3');
qp4ind     = findcells(ring,'FamName','QP4');
qp31ind     = findcells(ring,'FamName','QP31');
qp41ind     = findcells(ring,'FamName','QP41');
%pbind=findcells(r_original,'PolynomB');
focval_qp1 = atgetfieldvalues(ring,qp1ind,'PolynomB',{1,2});
focval_qp2 = atgetfieldvalues(ring,qp2ind,'PolynomB',{1,2});
focval_qp3 = atgetfieldvalues(ring,qp3ind,'PolynomB',{1,2});
focval_qp4 = atgetfieldvalues(ring,qp4ind,'PolynomB',{1,2});
focval_qp31 = atgetfieldvalues(ring,qp31ind,'PolynomB',{1,2});
focval_qp41 = atgetfieldvalues(ring,qp41ind,'PolynomB',{1,2});

qall = [focval_qp1 focval_qp2 focval_qp3 focval_qp4 focval_qp31 focval_qp41];