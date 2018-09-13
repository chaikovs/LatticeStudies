function varargout = TDR_good_017_064_r56_04_sx_Dff_BPM
%************************************************************
% ThomX ring lattice for AT (Accelrator Toolbox). 
%   Based on the Alexandre Loulergue's
%   Beta lattice file:
%     TDR_good_0.17_0.64_r56_0.2_sx_Dff.str,
%   which has the *** VERSION ***
%   BETA-LNS v8.00 /10/07/12/ 2012/12/12 14:48:0                                    
%CLS  - hexagone - tune=3.17  1.64 - C=16.80m h=28 - corr                       
%
%
% DON'T change the Lattice element names, since
% the lattice element names are fixed in the 
% MML file: /machine/THOMX/StorageRing/updateatindex.m
%
%  Written by Jianfeng Zhang @ LAL, 22/02/2013.
%  Modified by Iryna Chaikovska @ LAL, 27/02/2015.
%  tune=3.17  1.64 - C=17.9867m h=30 (BPM in good locations + COR)
%************************************************************
 
 global FAMLIST THERING GLOBVAL
 
 
 GLOBVAL.E0 = 50e6; % Ring energy [eV]
 GLOBVAL.LatticeFile = mfilename;
 FAMLIST = cell(0);
 
 disp(['** Loading THOMX ring lattice ', mfilename]);
 
 L0 = 17.986715999999987;%17.986704000000000;%17.986692; % design circumference length [m]; 
                          % 6D COD of the ThomX machine is very sensitive
                          %  to the accuracy of the circumference !!!!!!!!
                          %L0 =   findspos(THERING,length(THERING)+1)
 C0 = 2.99792458e8; % speed of light [m/s]

 HarmNumber = 30; % RF harmonic number

 %% ======================
%  RF Cavity
% ========================
%              NAME   L     U[V]       f[Hz]          h        method
RF = rfcavity('RF', 0.0,  300e3, HarmNumber*C0/L0, HarmNumber, ...
               'CavityPass');

%% =======================
% drift
%%========================
%sextupole length
LSX = 0;%0.06;



SD0 = drift('SD0', 0.1000000E+00, 'DriftPass');

SD1S1B = drift('SD1S1B', 0.600000E-01, 'DriftPass');
SD1S1 = drift('SD1S1', 0.190000E+00-LSX/2, 'DriftPass');
SD1S = drift('SD1S', 0.7500000E-01-LSX/2, 'DriftPass');

SD3S2= drift('SD3S2', 0.1250000E+00-LSX/2, 'DriftPass');

SD2 = drift('SD2', 0.2100000E+00, 'DriftPass');
SD31 = drift('SD31', 0.7537500E+00, 'DriftPass');
% SD31h = drift('SD31h', 0.7537500E+00/2, 'DriftPass');
% SD32K = drift('SD32K', 0.3600000E+00, 'DriftPass');
% SD32 = drift('SD32', 0.1937500E+00, 'DriftPass');
SD5B = drift('SD5B', 0.7400000E+00, 'DriftPass');
SD5 = drift('SD5', 0.0600000E+00, 'DriftPass');

SD3S = drift('SD3S', 0.7500000E-01-LSX/2, 'DriftPass');
SD3S1 = drift('SD3S1', 0.1350000E+00-LSX/2, 'DriftPass');

SD2SB = drift('SD2SB', 0.4000000E-01, 'DriftPass');
SD2S = drift('SD2S', 0.1950000E+00, 'DriftPass');

%% =======================
% dipole
%========================
L = 0.27646;
theta = 0.785398;
thetae = 0.0;
thetas = 0.0;
K =0.0;
beta_gap=0.04;
tracy_gap = beta_gap*2*0.348;
fullgap = tracy_gap;
edge_effect1 = 1;
edge_effect2 = 1;

% BEND  =  rbend3('BEND', L, theta, thetae, thetas, K,fullgap,edge_effect1,edge_effect2, ...
%                'BndMPoleSymplecticNew4Pass');
BEND  =  rbend3('BEND', L, theta, thetae, thetas, K,fullgap,edge_effect1,edge_effect2, ...
               'BendLinearPass');
 
 
 
%  COE = drift('COE',0.0,'DriftPass');
%  COS = drift('COS',0.0,'DriftPass');

%% =======================
% quadrupole
%========================
LQP = 0.15; % quadrupole length
QPassMethod = 'StrMPoleSymplectic4Pass'; % tracking method

QP1 = quadrupole('QP1', LQP, -.3044637E+01, QPassMethod);
QP2 = quadrupole('QP2', LQP, 0.7621329E+01, QPassMethod);
QP3 = quadrupole('QP3', LQP, -.1692568E+02, QPassMethod);
QP4 = quadrupole('QP4', LQP, 0.1806677E+02, QPassMethod);
QP31 = quadrupole('QP31', LQP, -.1013750E+02, QPassMethod);
QP41 = quadrupole('QP41', LQP, 0.2401124E+01, QPassMethod);

  
%% =======================
% sextupole
%========================

sx_on = 1; 

 SPassMethod = 'StrMPoleSymplectic4Pass';

%  SX1 = sextupole('SX1', LSX, -.6204968E+01/LSX*sx_on,  SPassMethod); 
%  SX2 = sextupole('SX2', LSX,  0.2381947E+01/LSX*sx_on, SPassMethod); 
%  SX3 = sextupole('SX3', LSX, -0.3115290E+01/LSX*sx_on, SPassMethod); 

 SX1 = sextupole('SX1', 0.1000000E-05, -.6204968E+07*sx_on,  SPassMethod); 
 SX2 = sextupole('SX2', 0.1000000E-05,  0.2381947E+07*sx_on, SPassMethod); 
 SX3 = sextupole('SX3', 0.1000000E-05, -0.3115290E+07*sx_on, SPassMethod); 

%% =======================
% BPM
% installed before or after quadrupole
%BPM = marker('BPM', 'IdentityPass');
%========================
% horizontal
BPMx = marker('BPMx', 'IdentityPass');
% vertical
BPMz = marker('BPMz', 'IdentityPass');
 
%% =======================
% correctors
% the same location as sextupole
%========================
% length, kick angle x/y
 HCOR = corrector('HCOR',1.0e-6,[0.0, 0],'CorrectorPass');
 
 VCOR = corrector('VCOR',1.0e-6,[0, 0.0],'CorrectorPass');

 %%==========================
 %   MARKER
 %%==========================
 
DEBUT =  marker('DEBUT','IdentityPass');
FIN   =  marker('FIN','IdentityPass');

%{**************}
%{  injection   }
%{**************}
 SEPT = marker('SEPT','IdentityPass');
 KICKER = marker('KICKer','IdentityPass');

%% Lattice definition
%=======================
% build lattice
%========================

str1 = [SD5B BPMx BPMz SD5 QP1 SD2 QP2 SD31];
cell1 = [ BEND SD3S2 SX1 HCOR VCOR SD3S QP31 SD3S SX3 HCOR VCOR SD3S1 QP41 SD1S1B BPMx BPMz SD1S1 SX2 HCOR VCOR SD1S QP4 SD2 QP3 SD2SB BPMx BPMz SD2S BEND SD0];
cell2 = [SD0 BEND SD2S BPMx BPMz SD2SB QP3 SD2 QP4 SD1S SX2 HCOR VCOR SD1S1 BPMx BPMz SD1S1B QP41 SD3S1 SX3 HCOR VCOR SD3S QP31 SD3S SX1 HCOR VCOR SD3S2 BEND];
str2 = [SD31 QP2 SD2 QP1 SD5 BPMx BPMz SD5B];


ELIST =[DEBUT RF str1 cell1 cell2 str2 SEPT str1 cell1 cell2 str2 FIN];
      
              
buildlat(ELIST);

% Set all magnets to same energy
THERING = setcellstruct(THERING,'Energy',1:length(THERING),GLOBVAL.E0);
    
evalin('caller','global THERING FAMLIST GLOBVAL');

% print the summary of the lattice
%atsummary;

if nargout
    varargout{1} = THERING;
end











  
