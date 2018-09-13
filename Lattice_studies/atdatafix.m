function [RING] = atdatafix(RING)
% Higher order chromaticities

RING=fittune2_alex(RING, [2.760 0.910], 'QF66E', 'B61HE');
RING=fitchrom_alex(RING,[-0.1 -0.12],'SXD2E' ,'SXF3E' );

end

