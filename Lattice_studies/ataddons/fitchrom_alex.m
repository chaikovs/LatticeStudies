function THERING =  fitchrom_alex(THERING,newchrom, sextfam1, sextfam2, varargin)
%FITCHROM2 fits chromaticity  of THERING using 2 sextupole families
%  FITCHROM2(NEWCHROM,SEXTUPOLEFAMILY1,SEXTUPOLEFAMILY2)
% 
%  INPUTS
%  1. newchrom - 2D chromaticity vector to fit to
%  2. sextfam1 - Family name for the first sextupole family
%  3. sextfam1 - Family name for the second sextupole family
%  4. Display  - Displays fitting results {default}
%     NoDisplay- Do not displays fitting results
% 
%  OUTPUTS
%  1. changes in sextupoles strength
%
%  EXAMPLES
%  1. fitchrom2([2 2],'S9','S10')
%
%  See Also fittune2

%
%  Written by Andrei Terebilo
%  Modified by Laurent S. Nadolski
%  MARCH 25, 2005 - Take into account thin sextupoles
%                 - Display Flag

DisplayFlag = 1;

%% Optional Input data parser
for i = length(varargin):-1:1
    if strcmpi(varargin{i},'Display')
        DisplayFlag = 1;
        varargin(i) = [];
    elseif strcmpi(varargin{i},'NoDisplay')
        DisplayFlag = 0;
        varargin(i) = [];
    end
end

% Must declare THERING as global in order for the function to modify sextupole values 
% if ~isempty(whos('global','THERING'))
%     global THERING
% end

%make a column vector
newchrom = newchrom(:);
% Thick sextupole
deltaS = 1e-3; % step size in Sextupole strength
deltaP = 1e-8;

% find indexes of the 2 sextupole families use for fitting
S1I = findcells(THERING,'FamName',sextfam1);
S2I = findcells(THERING,'FamName',sextfam2);
InitialS1 = getcellstruct(THERING,'PolynomB',S1I,3);
InitialS2 = getcellstruct(THERING,'PolynomB',S2I,3);

% Thin sextupoles
if THERING{S1I(1)}.Length < 1e-3 
  deltaS = 1e-2*1e8; % step size in Sextupole strength
end

% Compute initial tunes and chromaticities before fitting 

[ LD, InitialTunes] = linopt(THERING,0);
[ LDdP, ITdP] =linopt(THERING,deltaP);

InitialChrom = (ITdP-InitialTunes)/deltaP;

TempTunes = InitialTunes;
TempChrom = InitialChrom;
TempS1 = InitialS1; 
TempS2 = InitialS2;

%% loop over n times to assure convergence
for i=1:5
		
	% Take Derivative
	THERING = setcellstruct(THERING,'PolynomB',S1I,TempS1+deltaS,3);
	[LD , Tunes_dS1 ] = linopt(THERING,0);
	[LD , Tunes_dS1dP ] = linopt(THERING,deltaP);

	THERING = setcellstruct(THERING,'PolynomB',S1I,TempS1,3);
	THERING = setcellstruct(THERING,'PolynomB',S2I,TempS2+deltaS,3);
	[LD , Tunes_dS2 ] = linopt(THERING,0);
	[LD , Tunes_dS2dP ] = linopt(THERING,deltaP);
	THERING = setcellstruct(THERING,'PolynomB',S2I,TempS2,3);

	%Construct the Jacobian
	Chrom_dS1 = (Tunes_dS1dP-Tunes_dS1)/deltaP;
	Chrom_dS2 = (Tunes_dS2dP-Tunes_dS2)/deltaP;

	J = ([Chrom_dS1(:) Chrom_dS2(:)] - [TempChrom(:) TempChrom(:)])/deltaS;
	Jinv = inv(J);

	dchrom = (newchrom(:) - TempChrom(:));
	dS = Jinv*dchrom;

	TempS1 = TempS1 + dS(1);
	TempS2 = TempS2 + dS(2);

	THERING = setcellstruct(THERING,'PolynomB',S1I,TempS1,3);
	THERING = setcellstruct(THERING,'PolynomB',S2I,TempS2,3);

	[ LD, TempTunes] = linopt(THERING,0);
	[ LD, TempTunesdP] = linopt(THERING,deltaP);
	TempChrom = (TempTunesdP-TempTunes)/deltaP;

end

% %% Display how good is the fit
% if DisplayFlag
%     fprintf('Desired chromaticities xix=%f xiz=%f\n',newchrom);
%     [tune xi] = tunechrom(THERING,0,'chrom');
%     fprintf('Reached chromaticities xix=%f xiz=%f\n',xi);
% end