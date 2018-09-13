function z=rbend3(fname,L,A,A1,A2,K,gap,EdgeEffect1,EdgeEffect2, method)
%RBEND2('FAMILYNAME',  Length[m], BendingAngle[rad], EntranceAngle[rad],
%	ExitAngle[rad], K, gap, F1, F2, 'METHOD')
%	creates a new family in the FAMLIST - a structure with fields
%     FamName        	family name
%     Length         	length of the arc for an on-energy particle [m]
%     BendingAngle	total bending angle [rad]
%     EntranceAngle	[rad] (0 - for sector bends)
%     ExitAngle		[rad] (0 - for sector bends)
%     ByError	        error in the dipole field relative to the design value 
%     K			quadrupole K-value for combined funtion bends
%     gap               FullGap
%     EdgeEffect1       flag to turn on (1) or off (0) the dipole
%                       fringe field and edge focus at the entrance
%                       of the dipole
%     EdgeEffect2       flag to turn on (1) or off (0) the dipole
%                       fringe field and edge focus at the exit
%                       of the dipole
%     PassMethod        name of the function to use for tracking
% returns assigned address in the FAMLIST that is uniquely identifies
% the family
% 
% Added by Laurent S. Nadolski, SOLEIL, 03/04
%
%
% Modified by Jianfeng Zhang @ LAL, 29/04/2013
%  to fix the bug to turn on/off dipole FF and edge effects in
%  sector dipoles.

ElemData.Class      	= 'Bend';
ElemData.FamName = fname;  % add check for identical family names
ElemData.Length			= L;
ElemData.MaxOrder		= 3;
ElemData.NumIntSteps 	= 10;
ElemData.BendingAngle  	= A;
ElemData.EntranceAngle 	= A1;
ElemData.ExitAngle     	= A2;
ElemData.ByError     	= 0;
ElemData.K      		= K;
ElemData.FullGap   		= gap;
ElemData.FringeInt1	    = 0.5; % same convention as in Tracy II
ElemData.FringeInt2	    = 0.5; % same convention as in Tracy II
ElemData.EdgeEffect1    = EdgeEffect1;% flag to turn on (1) or off (0) edge 
                                      % effects and fringe field at the 
                                      % entrance of the dipole. 
                                      % Added by Jianfeng Zhang @ LAL,
                                      % 29/04/2013
ElemData.EdgeEffect2    = EdgeEffect2;% flag to turn on (1) or off (0) edge 
                                      % effects and fringe field at the
                                      % exit of the dipole.
                                      % Added by Jianfeng Zhang @ LAL,
                                      % 29/04/2013


ElemData.R1             = diag(ones(6,1)); %rotation error
ElemData.R2             = diag(ones(6,1));
ElemData.T1             = zeros(1,6); %displacement error
ElemData.T2             = zeros(1,6);

ElemData.PolynomA		= [0 0 0 0];	%skew component of field 
ElemData.PolynomB		= [0 K 0 0];    %normal component of field
ElemData.PassMethod 	= method;

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName  = fname;
FAMLIST{z}.NumKids  = 0;
FAMLIST{z}.KidsList = [];
FAMLIST{z}.ElemData = ElemData;
