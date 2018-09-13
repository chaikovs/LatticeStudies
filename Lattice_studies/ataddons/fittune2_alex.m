function THERING =  fittune2_alex(THERING, newtunes, quadfam1, quadfam2, varargin)
%FITTUNE2 fits linear tunes of THERING using 2 quadrupole families
% FITTUNE2(NEWTUNES,QUADFAMILY1,QUADFAMILY2)
%  INPUTS
%  1. newtunes - 2D tune vector to fit to
%  2. quadfam1 - Family name for the first quadrupole family
%  3. quadfam1 - Family name for the second quadrupole family
%  4. delta    - Kvariation for computing Jacobian matrix 
%  5. Display  - Displays fitting results {default}
%     NoDisplay- Do not displays fitting results
%
%  EXAMPLES
%  1. fittune2([0.2 0.3],'Q7','Q9')
%
%  See Also fitchrom2

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

if length(varargin) > 3 % use externally supplied step size for quadrupole K-values 
    delta = varargin{1}
else
    delta = 1e-6; % default step size for quadrupole K-values 
end
% find indexes of the 2 quadrupole families use for fitting
Q1I = findcells(THERING,'FamName',quadfam1);
Q2I = findcells(THERING,'FamName',quadfam2);

%InitialK1 = getcellstruct(THERING,'K',Q1I);
%InitialK2 = getcellstruct(THERING,'K',Q2I);
InitialPolB1 = getcellstruct(THERING,'PolynomB',Q1I,2);
InitialPolB2 = getcellstruct(THERING,'PolynomB',Q2I,2);

%% Compute initial tunes before fitting 
%[ LD, InitialTunes] = linopt(THERING,0);
[~,  InitialTunes] = twissring(THERING, 0, 1:length(THERING)+1);

TempTunes = InitialTunes;
%TempK1 = InitialK1;
%TempK2 = InitialK2;
TempPolB1 = InitialPolB1;
TempPolB2 = InitialPolB2;

%% Take Derivative
%THERING = setcellstruct(THERING,'K',Q1I,TempK1+delta);
THERING = setcellstruct(THERING,'PolynomB',Q1I,TempPolB1+delta,2);
[~ , Tunes_dK1 ] = twissring(THERING, 0, 1:length(THERING)+1);
%THERING = setcellstruct(THERING,'K',Q1I,TempK1);
THERING = setcellstruct(THERING,'PolynomB',Q1I,TempPolB1,2);
%THERING = setcellstruct(THERING,'K',Q2I,TempK2+delta);
THERING = setcellstruct(THERING,'PolynomB',Q2I,TempPolB2+delta,2);
[~ , Tunes_dK2 ] = twissring(THERING, 0, 1:length(THERING)+1);
%THERING = setcellstruct(THERING,'K',Q2I,TempK2);
THERING = setcellstruct(THERING,'PolynomB',Q2I,TempPolB2,2);

%% Construct the Jacobian
J = ([Tunes_dK1(:) Tunes_dK2(:)] - [TempTunes(:) TempTunes(:)])/delta;
Jinv = inv(J);

dnu = (newtunes(:) - TempTunes(:));
dK = Jinv*dnu;

%TempK1 = TempK1+dK(1);
%TempK2 = TempK2+dK(2);
TempPolB1 = TempPolB1 + dK(1);
TempPolB2 = TempPolB2 + dK(2);

%THERING = setcellstruct(THERING,'K',Q1I,TempK1);
THERING = setcellstruct(THERING,'PolynomB',Q1I,TempPolB1,2);
%THERING = setcellstruct(THERING,'K',Q2I,TempK2);
THERING = setcellstruct(THERING,'PolynomB',Q2I,TempPolB2,2);

%[LD,TempTunes] = linopt(THERING,0);
[~,TempTunes ] = twissring(THERING, 0, 1:length(THERING)+1);

%InitialK1 - TempK1;
%InitialK2 - TempK2;
%TempTunes

%% Display how good is the fit
if DisplayFlag
    fprintf('Desired tunes nux=%f nuz=%f\n',newtunes);
    [tune xi] = tunechrom(THERING,0,'chrom');
    fprintf('Reached tunes nux=%f nuz=%f\n',TempTunes);
end