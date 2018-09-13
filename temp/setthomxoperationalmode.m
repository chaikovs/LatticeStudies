function setthomxoperationalmode(ModeNumber)
%SETOPERATIONALMODE - Switches between the various operational modes
%  setoperationalmode(ModeNumber)
%
%  INPUTS
%  1. ModeNumber = number
%        1 '50 MeV, 3.17 1.72, r56=0.2', ...

%  See also setoperationalmode, aoinit, updateatindex, soleilinit, setmmldirectories, lattice_prep

%  NOTES
% use local_set_config_mode for defining status of S11 et S12;


% Based on setoperationalmode.m, 
% Modified by Jianfeng Zhang @ LAL 18/06/2013
% FOR CVS
% $Header$
global THERING

% Check if the AO exists
checkforao;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerator Dependent Modes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    ModeNumber = [];
end
if isempty(ModeNumber)
    ModeCell = {...
        ' 1/ 50 MeV, 3.17 1.64, r56=0.2', ...
        };
    [ModeNumber, OKFlag] = listdlg('Name','THOMX','PromptString','Select the Operational Mode:', ...
        'SelectionMode','single', 'ListString', ModeCell, 'ListSize', [450 200], 'InitialValue', 1);
    if OKFlag ~= 1
        fprintf('   Operational mode not changed\n');
        return
   % elseif ModeNumber == length(ModeCell);
   %     ModeNumber = 100;  % Laurent
   % end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accelerator Data Structure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AD = getad;
AD.Machine = 'THOMX';            % Will already be defined if setpathmml was used
AD.MachineType = 'StorageRing';   % Will already be defined if setpathmml was used
AD.SubMachine  = 'StorageRing';   % Will already be defined if setpathmml was used
AD.OperationalMode = '';          % Gets filled in later
AD.HarmonicNumber = 30;

% Defaults RF for dispersion and chromaticity measurements (must be in Hardware units)
AD.DeltaRFDisp = 100e-6;
AD.DeltaRFChro = [-100 -50 0 50 100] * 1e-6;

% Tune processor delay: delay required to wait
% to have a fresh tune measurement after changing
% a variable like the RF frequency.  Setpv will wait
% 2.2 * TuneDelay to be guaranteed a fresh data point.
%AD.BPMDelay  = 0.25; % use [N, BPMDelay]=getbpmaverages (AD.BPMDelay will disappear)
AD.TuneDelay = 1;

% The offset and golden orbits are stored at the end of this file
% TODO
%BuildOffsetAndGoldenOrbits;  % Local function


% SP-AM Error level
% AD.ErrorWarningLevel = 0 -> SP-AM errors are Matlab errors {Default}
%                       -1 -> SP-AM errors are Matlab warnings
%                       -2 -> SP-AM errors prompt a dialog box
%                       -3 -> SP-AM errors are ignored (ErrorFlag=-1 is returned)
AD.ErrorWarningLevel = 0;

%%%%%%%%%%%%%%%%%%%%%
% Operational Modes %
%%%%%%%%%%%%%%%%%%%%%

% Mode setup variables (mostly path and file names)
% AD.OperationalMode - String used in titles
% ModeName - String used for mode directory name off DataRoot/MachineName
% OpsFileExtension - string add to default file names

 


%% ModeNumber == 9 bx=10m nominal lattice 2010 until installation of S11 and S12
if ModeNumber == 9
    % User mode - nominal lattice 2010 until installation of S11 and S12
    AD.OperationalMode = '2.7391 GeV, 18.2 10.3';
    AD.Energy = 2.7391; % Make sure this is the same as bend2gev at the production lattice!
    ModeName = 'lat_2020_3170b';
    OpsFileExtension = '_lat_2020_3170b';

    % AT lattice
    AD.ATModel = 'lat_2020_3170b';
    eval(AD.ATModel);  %run model for compiler;

    % Golden TUNE is with the TUNE family
    % 18.2020 / 10.3170
    AO = getao;
    AO.TUNE.Monitor.Golden = [
        0.2020
        0.3170
        NaN];

    % Golden chromaticity is in the AD (Physics units)
    AD.Chromaticity.Golden = [2; 2];

    % Status factory
    local_set_config_mode('normalconfig120');
    AO = local_setmagnetcoefficient(AO, @magnetcoefficients_new_calib_new_modele_juin2009);
    setao(AO);

    %% ModeNumber == 19 User mode - S11 betax=10m till November 2010
elseif ModeNumber == 19 
    % User mode - S11 betax=10m 2010
    AD.OperationalMode = '2.7391 GeV, 18.202 10.317 S11';
    AD.Energy = 2.7391; % Make sure this is the same as bend2gev at the production lattice!
    ModeName = 'lat_2020_3170f';
    OpsFileExtension = '_lat_2020_3170f';

    % AT lattice
    AD.ATModel = 'lat_2020_3170f';
    eval(AD.ATModel);  %run model for compiler;

    % Golden TUNE is with the TUNE family
    % 18.2020 / 10.3170
    AO = getao;
    AO.TUNE.Monitor.Golden = [
        0.2020
        0.3170
        NaN];

    % Golden chromaticity is in the AD (Physics units)
    AD.Chromaticity.Golden = [2; 2.6];

    % Status factory
    local_set_config_mode('S11config120');
    AO = local_setmagnetcoefficient(AO, @magnetcoefficients_new_calib_new_modele_juin2009);
    setao(AO);
   
    %% ModeNumber == 16 User mode - betax = 5m
elseif ModeNumber == 16
    % User mode - betax = 5m
    AD.OperationalMode = '2.7391 GeV, 18.2 10.3';
    AD.Energy = 2.7391; % Make sure this is the same as bend2gev at the production lattice!
    ModeName = 'lat_2020_3170e';
    OpsFileExtension = '_lat_2020_3170e';

    % AT lattice
    AD.ATModel = 'lat_2020_3170e';
    eval(AD.ATModel);  %run model for compiler;

    % Golden TUNE is with the TUNE family
    % 18.2020 / 10.3170
    % 18.1990 / 10.3170 April 2011
    AO = getao;
    AO.TUNE.Monitor.Golden = [
        0.1990
        0.3100
        0.00642];

    % Golden chromaticity is in the AD (Physics units)
    AD.Chromaticity.Golden = [2; 2];

    local_set_config_mode('S11config120');
    AO = local_setmagnetcoefficient(AO, @magnetcoefficients_new_calib_new_modele_juin2009);
    setao(AO);

    %% ModeNumber == 10 User mode - with PX2 corrector   
elseif ModeNumber == 10
    % User mode - with PX2 corrector
    AD.OperationalMode = '2.7391 GeV, 18.2 10.3';
    AD.Energy = 2.7391; % Make sure this is the same as bend2gev at the production lattice!
    ModeName = 'lat_2020_3170bPX2';
    OpsFileExtension = '_lat_2020_3170bPX2';

    % AT lattice
    AD.ATModel = 'lat_2020_3170bPX2';
    eval(AD.ATModel);  %run model for compiler;

    % Golden TUNE is with the TUNE family
    % 18.2020 / 10.3170
    AO = getao;
    AO.TUNE.Monitor.Golden = [
        0.2020
        0.3170
        NaN];

    % Golden chromaticity is in the AD (Physics units)
    AD.Chromaticity.Golden = [2; 2];

    local_set_config_mode('normalconfig120');
    AO = local_setmagnetcoefficient(AO, @magnetcoefficients_new_calib_new_modele_juin2009);
    setao(AO);

    %% ModeNumber == 11  User mode - Nanoscopium
elseif ModeNumber == 11 % User mode - Nanoscopium
    AD.OperationalMode = '2.7391 GeV, 18.2175 10.3120';
    AD.Energy = 2.7391; % Make sure this is the same as bend2gev at the production lattice!
    ModeName = 'nano_2175_3120a';
    OpsFileExtension = '_nano_2175_3120a';

    % AT lattice
    AD.ATModel = 'nano_2175_3120'; % new lattice version from Alex
    eval(AD.ATModel);  %run model for compiler;

    % Golden TUNE is with the TUNE family
    AO = getao;
    AO.TUNE.Monitor.Golden = [
        0.2175
        0.3120
        NaN];

    % Golden chromaticity is in the AD (Physics units)
    AD.Chromaticity.Golden = [2; 2];

    local_set_config_mode('nanoscopiumconfig');
    AO = local_setmagnetcoefficient(AO, @magnetcoefficients_new_calib_new_modele_juin2009);
    setao(AO);

  
    % triplet upstreams and downstreams of SDL13 for nanoscopium
    % Need to point to another family for magnetcoefficients (other range of current)
    % Q1 upstream
    HW2PhysicsParams = magnetcoefficients_new_calib_new_modele_juin2009('Q1');
    AO.Q1.Monitor.HW2PhysicsParams{1}(6,:)  = HW2PhysicsParams;
    AO.Q1.Monitor.Physics2HWParams{1}(6,:)  = HW2PhysicsParams;
    AO.Q1.Setpoint.HW2PhysicsParams{1}(6,:) = HW2PhysicsParams;
    AO.Q1.Setpoint.Physics2HWParams{1}(6,:)  = HW2PhysicsParams;
    % Q1 downstream
    AO.Q1.Monitor.HW2PhysicsParams{1}(7,:)  = HW2PhysicsParams;
    AO.Q1.Monitor.Physics2HWParams{1}(7,:)  = HW2PhysicsParams;
    AO.Q1.Setpoint.HW2PhysicsParams{1}(7,:) = HW2PhysicsParams;
    AO.Q1.Setpoint.Physics2HWParams{1}(7,:) = HW2PhysicsParams;
    % Q2 upstream
    HW2PhysicsParams = magnetcoefficients_new_calib_new_modele_juin2009('Q2');
    AO.Q2.Monitor.HW2PhysicsParams{1}(6,:)  = HW2PhysicsParams;
    AO.Q2.Monitor.Physics2HWParams{1}(6,:)  = HW2PhysicsParams;
    AO.Q2.Setpoint.HW2PhysicsParams{1}(6,:) = HW2PhysicsParams;
    AO.Q2.Setpoint.Physics2HWParams{1}(6,:) = HW2PhysicsParams;
    % Q2 downstream
    AO.Q2.Monitor.HW2PhysicsParams{1}(7,:)  = HW2PhysicsParams;
    AO.Q2.Monitor.Physics2HWParams{1}(7,:)  = HW2PhysicsParams;
    AO.Q2.Setpoint.HW2PhysicsParams{1}(7,:) = HW2PhysicsParams;
    AO.Q2.Setpoint.Physics2HWParams{1}(7,:) = HW2PhysicsParams;
    % Q3 upstream
    HW2PhysicsParams = magnetcoefficients_new_calib_new_modele_juin2009('Q1');
    AO.Q3.Monitor.HW2PhysicsParams{1}(6,:)  = HW2PhysicsParams;
    AO.Q3.Monitor.Physics2HWParams{1}(6,:)  = HW2PhysicsParams;
    AO.Q3.Setpoint.HW2PhysicsParams{1}(6,:) = HW2PhysicsParams;
    AO.Q3.Setpoint.Physics2HWParams{1}(6,:) = HW2PhysicsParams;
    % Q3 downstream
    AO.Q3.Monitor.HW2PhysicsParams{1}(7,:)  = HW2PhysicsParams;
    AO.Q3.Monitor.Physics2HWParams{1}(7,:)  = HW2PhysicsParams;
    AO.Q3.Setpoint.HW2PhysicsParams{1}(7,:) = HW2PhysicsParams;
    AO.Q3.Setpoint.Physics2HWParams{1}(7,:) = HW2PhysicsParams;
    
    setao(AO);

    %% ModeNumber == 14 % User mode - Nanoscopium
elseif ModeNumber == 14 % User mode - Nanoscopium
    AD.OperationalMode = '2.7391 GeV, 18.2175 10.3120';
    AD.Energy = 2.7391; % Make sure this is the same as bend2gev at the production lattice!
    ModeName = 'nano_8000_7000a';
    OpsFileExtension = '_nano_8000_7000a';

    % AT lattice
    AD.ATModel = 'nano_8000_7000'; % new lattice version from Alex
    eval(AD.ATModel);  %run model for compiler;

    % Golden TUNE is with the TUNE family
    AO = getao;
    AO.TUNE.Monitor.Golden = [
        0.2175
        0.3120
        NaN];

    % Golden chromaticity is in the AD (Physics units)
    AD.Chromaticity.Golden = [2; 2];


    local_set_config_mode('nanoscopiumconfig');
    AO = local_setmagnetcoefficient(AO, @magnetcoefficients_new_calib_new_modele_juin2009);
    setao(AO);

    % triplet upstreams and downstreams of SDL13 for nanoscopium
    % Need to point to another family for magnetcoefficients (other range of current)
    % Q1 upstream
    HW2PhysicsParams = magnetcoefficients_new_calib_new_modele_juin2009_nano_80_70('Q1');
    AO.Q1.Monitor.HW2PhysicsParams{1}(6,:)  = HW2PhysicsParams;
    AO.Q1.Monitor.Physics2HWParams{1}(6,:)  = HW2PhysicsParams;
    AO.Q1.Setpoint.HW2PhysicsParams{1}(6,:) = HW2PhysicsParams;
    AO.Q1.Setpoint.Physics2HWParams{1}(6,:)  = HW2PhysicsParams;
    % Q1 downstream
    AO.Q1.Monitor.HW2PhysicsParams{1}(7,:)  = HW2PhysicsParams;
    AO.Q1.Monitor.Physics2HWParams{1}(7,:)  = HW2PhysicsParams;
    AO.Q1.Setpoint.HW2PhysicsParams{1}(7,:) = HW2PhysicsParams;
    AO.Q1.Setpoint.Physics2HWParams{1}(7,:) = HW2PhysicsParams;
    % Q2 upstream
    HW2PhysicsParams = magnetcoefficients_new_calib_new_modele_juin2009_nano_80_70('Q2');
    AO.Q2.Monitor.HW2PhysicsParams{1}(6,:)  = HW2PhysicsParams;
    AO.Q2.Monitor.Physics2HWParams{1}(6,:)  = HW2PhysicsParams;
    AO.Q2.Setpoint.HW2PhysicsParams{1}(6,:) = HW2PhysicsParams;
    AO.Q2.Setpoint.Physics2HWParams{1}(6,:) = HW2PhysicsParams;
    % Q2 downstream
    AO.Q2.Monitor.HW2PhysicsParams{1}(7,:)  = HW2PhysicsParams;
    AO.Q2.Monitor.Physics2HWParams{1}(7,:)  = HW2PhysicsParams;
    AO.Q2.Setpoint.HW2PhysicsParams{1}(7,:) = HW2PhysicsParams;
    AO.Q2.Setpoint.Physics2HWParams{1}(7,:) = HW2PhysicsParams;
    % Q3 upstream
    HW2PhysicsParams = magnetcoefficients_new_calib_new_modele_juin2009_nano_80_70('Q1');
    AO.Q3.Monitor.HW2PhysicsParams{1}(6,:)  = HW2PhysicsParams;
    AO.Q3.Monitor.Physics2HWParams{1}(6,:)  = HW2PhysicsParams;
    AO.Q3.Setpoint.HW2PhysicsParams{1}(6,:) = HW2PhysicsParams;
    AO.Q3.Setpoint.Physics2HWParams{1}(6,:) = HW2PhysicsParams;
    % Q3 downstream
    AO.Q3.Monitor.HW2PhysicsParams{1}(7,:)  = HW2PhysicsParams;
    AO.Q3.Monitor.Physics2HWParams{1}(7,:)  = HW2PhysicsParams;
    AO.Q3.Setpoint.HW2PhysicsParams{1}(7,:) = HW2PhysicsParams;
    AO.Q3.Setpoint.Physics2HWParams{1}(7,:) = HW2PhysicsParams;

    % triplet nanoscopium
    HW2PhysicsParams = magnetcoefficients_new_calib_new_modele_juin2009_nano_80_70('Q11');
    AO.Q11.Monitor.HW2PhysicsParams{1}(1,:)  = HW2PhysicsParams;
    AO.Q11.Setpoint.HW2PhysicsParams{1}(1,:) = HW2PhysicsParams;
    AO.Q11.Monitor.Physics2HWParams{1}(1,:)  = HW2PhysicsParams;
    AO.Q11.Setpoint.Physics2HWParams{1}(1,:) = HW2PhysicsParams;
    AO.Q11.Monitor.HW2PhysicsParams{1}(2,:)  = HW2PhysicsParams;
    AO.Q11.Setpoint.HW2PhysicsParams{1}(2,:) = HW2PhysicsParams;
    AO.Q11.Monitor.Physics2HWParams{1}(2,:)  = HW2PhysicsParams;
    AO.Q11.Setpoint.Physics2HWParams{1}(2,:) = HW2PhysicsParams;

    HW2PhysicsParams = magnetcoefficients_new_calib_new_modele_juin2009_nano_80_70('Q2');
    AO.Q12.Monitor.HW2PhysicsParams{1}(1,:)  = HW2PhysicsParams;
    AO.Q12.Setpoint.HW2PhysicsParams{1}(1,:) = HW2PhysicsParams;
    AO.Q12.Monitor.Physics2HWParams{1}(1,:)  = HW2PhysicsParams;
    AO.Q12.Setpoint.Physics2HWParams{1}(1,:) = HW2PhysicsParams;


    AO.Q11.Status = [1; 1];
    AO.Q12.Status = 1;

    setao(AO);
    
    %% ModeNumber == 1
elseif ModeNumber == 1
    % User mode - 
    AD.OperationalMode = '2.7391 GeV, 18.2 10.3';
    AD.Energy = 2.7391; % Make sure this is the same as bend2gev at the production lattice!
    ModeName = 'solamor2c';
    OpsFileExtension = '_solamor2c';

    % AT lattice
    AD.ATModel = 'solamor2linc';
    eval(AD.ATModel);  %run model for compilersolamor2linb;

    % Golden TUNE is with the TUNE family
    % 18.2020 / 10.3170
    AO = getao;
    AO.TUNE.Monitor.Golden = [
        0.2020
        0.3170
        NaN];
    % Golden chromaticity is in the AD (Physics units)
    AD.Chromaticity.Golden = [2; 2];

    % Status factory
    local_set_config_mode('normalconfig120');
    AO = local_setmagnetcoefficient(AO, @magnetcoefficients);
    setao(AO);

end
    
    
% Force units to hardware
switch2hw;

% Activation of correctors of HU640
if ModeNumber == 6
  switchHU640Cor('ON');
else
  switchHU640Cor('OFF');
end

% Set the AD directory path
setad(AD);
MMLROOT = setmmldirectories(AD.Machine, AD.SubMachine, ModeName, OpsFileExtension);
AD = getad;


% SOLEIL specific path changes

% Top Level Directories

%AD.Directory.DataRoot       = fullfile(MMLROOT, 'measdata', 'SOLEIL', 'StorageRingdata', filesep);
% RUCHE
MMLDATAROOT = getmmldataroot;
AD.Directory.Lattice        = fullfile(MMLROOT, 'machine', 'SOLEIL', 'StorageRing', 'Lattices', filesep);
AD.Directory.Orbit          = fullfile(MMLROOT, 'machine', 'SOLEIL', 'StorageRing',  'orbit', filesep);

% Data Archive Directories DO NOT REMOVE LINES

% Insertion Devices


% STANDALONE matlab applications
AD.Directory.Standalone     = fullfile(MMLROOT, 'machine', 'THOMX', 'standalone_applications', filesep);

% FOFB matlab applications

% For coupling correction. Used by coupling.m

% AD.Directory.InterlockData  = fullfile(AD.Directory.DataRoot, 'Interlock/'];

%Response Matrix Directories

% used by energytunette

% used by MAT's Steerette application


% Postmortem DATA

%Default Data File Prefix
AD.Default.BPMArchiveFile      = 'BPM';                %file in AD.Directory.BPM               orbit data
AD.Default.TuneArchiveFile     = 'Tune';               %file in AD.Directory.Tune              tune data
AD.Default.ChroArchiveFile     = 'Chro';               %file in AD.Directory.Chromaticity       chromaticity data
AD.Default.DispArchiveFile     = 'Disp';               %file in AD.Directory.Dispersion       dispersion data
AD.Default.CNFArchiveFile      = 'CNF';                %file in AD.Directory.CNF               configuration data
AD.Default.QUADArchiveFile     = 'QuadBeta';           %file in AD.Directory.QUAD             betafunction for quadrupoles
AD.Default.PINHOLEArchiveFile  = 'Pinhole';            %file in AD.Directory.PINHOLE             pinhole data
AD.Default.SkewArchiveFile     = 'SkewQuad';           %file in AD.Directory.SkewQuad             SkewQuad data
AD.Default.BBAArchiveFile      = 'BBA_DKmode';         %file in AD.Directory.BBA             BBA DK mode data

%Default Response Matrix File Prefix
AD.Default.BPMRespFile      = 'BPMRespMat';         %file in AD.Directory.BPMResponse       BPM response matrices
AD.Default.TuneRespFile     = 'TuneRespMat';        %file in AD.Directory.TuneResponse      tune response matrices
AD.Default.ChroRespFile     = 'ChroRespMat';        %file in AD.Directory.ChroResponse      chromaticity response matrices
AD.Default.DispRespFile     = 'DispRespMat';        %file in AD.Directory.DispResponse      dispersion response matrices
AD.Default.SkewRespFile     = 'SkewRespMat';        %file in AD.Directory.SkewResponse      skew quadrupole response matrices

%Orbit Control and Feedback Files
AD.Restore.GlobalFeedback   = 'Restore.m';

% Circumference
AD.Circumference = findspos(THERING,length(THERING)+1);
setad(AD);

% Updates the AT indices in the MiddleLayer with the present AT lattice
updateatindex;

% Set the model energy
setenergymodel(AD.Energy);


% Momentum compaction factor
MCF = getmcf('Model');
if isnan(MCF)
    AD.MCF = 4.498325442923014e-04;
    fprintf('   Model alpha calculation failed, middlelayer alpha set to  %f\n', AD.MCF);
else
    AD.MCF = MCF;
    fprintf('   Middlelayer alpha set to %f (AT model).\n', AD.MCF);
end
setad(AD);


% Add Gain & Offsets for magnet family
fprintf('   Setting magnet monitor gains based on the production lattice.\n');
%setgainsandoffsets;

%% Config texttalker (right location ?)
AD.TANGO.TEXTTALKERS={'ans/ca/texttalker.1', 'ans/ca/texttalker.2'};

% set LOCO gain and roll to zero
setlocodata('Nominal');

%%%%%%%%%%%%%%%%%%%%%%
% Final mode changes %
%%%%%%%%%%%%%%%%%%%%%%
if any(ModeNumber == [99])
    % User mode - 2.75 GeV, Nominal lattice

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add LOCO Parameters to AO and AT-Model %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     'Nominal'    - Sets nominal gains (1) / rolls (0) to the model.
    %     'SetGains'   - Set gains/coupling from a LOCO file.
    %     'SetModel'   - Set the model from a LOCO file.  But it only changes
    %                    the part of the model that does not get corrected
    %                    in 'Symmetrize' (also does a SetGains).
    %     'LOCO2Model' - Set the model from a LOCO file (also does a SetGains).
    %                    This sets all lattice machines fit in the LOCO run to the model.
    %
    % Basically, use 'SetGains' or 'SetModel' if the LOCO run was applied to the accelerator
    %            use 'LOCO2Model' if the LOCO run was made after the final setup

    % Store the LOCO file in the opsdata directory

    % MCF depends on optics !!!

    AD.OpsData.LOCOFile = [getfamilydata('Directory','OpsData'),'LOCO_163Quads_122BPMs'];
    
    try % TO BE DONE LATER IN 2012
        setlocodata('LOCO2Model', AD.OpsData.LOCOFile);
    catch
        fprintf('\n%s\n\n', lasterr);
        fprintf('   WARNING: there was a problem calibrating the model based on LOCO file %s.\n', AD.OpsData.LOCOFile);
    end

else
    setlocodata('Nominal');
end

fprintf('   lattice files have changed or if the AT lattice has changed.\n');
fprintf('   Middlelayer setup for operational mode: %s\n', AD.OperationalMode);

setad(orderfields(AD));

end

function local_setmagnetcoefficient(magnetcoeff_function)
% quadrupole magnet coefficients
% number of status 1 quadrupole families

AO = getao;
    
quadFamList = {'Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', ...
         'Q7', 'Q8', 'Q9', 'Q10'};
    
if family2status('Q11',1),
        quadFamList = [quadFamList, {'Q11'}];
end
    
if family2status('Q12',1),
        quadFamList = [quadFamList, {'Q12'}];
end
        
    
for k = 1:length(quadFamList),
        ifam = quadFamList{k};

        HW2PhysicsParams = feval(magnetcoeff_function, AO.(ifam).FamilyName);
        Physics2HWParams = HW2PhysicsParams;

        nb = size(AO.(ifam).DeviceName,1);

        for ii=1:nb,
            val = 1.0;
            AO.(ifam).Monitor.HW2PhysicsParams{1}(ii,:)                 = HW2PhysicsParams;
            AO.(ifam).Monitor.HW2PhysicsParams{2}(ii,:)                 = val;
            AO.(ifam).Monitor.Physics2HWParams{1}(ii,:)                 = Physics2HWParams;
            AO.(ifam).Monitor.Physics2HWParams{2}(ii,:)                 = val;
            AO.(ifam).Setpoint.HW2PhysicsParams{1}(ii,:)                = HW2PhysicsParams;
            AO.(ifam).Setpoint.HW2PhysicsParams{2}(ii,:)                = val;
            AO.(ifam).Setpoint.Physics2HWParams{1}(ii,:)                = Physics2HWParams;
            AO.(ifam).Setpoint.Physics2HWParams{2}(ii,:)                = val;
        end
end

% sextupole magnet coefficients
% number of status 1 sextupole families
sextuFamList = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6', ...
         'S7', 'S8', 'S9', 'S10'};
    
if family2status('S11',1),
        sextuFamList = [sextuFamList, {'S11'}];
end
    
if family2status('S12',1),
        sextuFamList = [sextuFamList, {'S12'}];
end

for k = 1:length(sextuFamList),
        ifam = sextuFamList{k};
        
        HW2PhysicsParams = feval(magnetcoeff_function, AO.(ifam).FamilyName);
        Physics2HWParams = HW2PhysicsParams;
        
        val = 1.0;
        AO.(ifam).Monitor.HW2PhysicsParams{1}(1,:)                 = HW2PhysicsParams;
        AO.(ifam).Monitor.HW2PhysicsParams{2}(1,:)                 = val;
        AO.(ifam).Monitor.Physics2HWParams{1}(1,:)                 = Physics2HWParams;
        AO.(ifam).Monitor.Physics2HWParams{2}(1,:)                 = val;
        AO.(ifam).Setpoint.HW2PhysicsParams{1}(1,:)                 = HW2PhysicsParams;
        AO.(ifam).Setpoint.HW2PhysicsParams{2}(1,:)                 = val;
        AO.(ifam).Setpoint.Physics2HWParams{1}(1,:)                 = Physics2HWParams;
        AO.(ifam).Setpoint.Physics2HWParams{2}(1,:)                 = val;
end
setao(AO);

end


function  local_set_config_mode(configmode)
% Function for activating new families of quadrupole and sextupoles
% magnets.

switch(configmode)
    case 'S11config120' % with S11 120 BPMs to be obsolete
        setfamilydata(1, 'S11', 'Status')
        setfamilydata(0, 'S12', 'Status')
        setfamilydata(0, 'Q11', 'Status')
        setfamilydata(0, 'Q12', 'Status')
        setfamilydata(0, 'HCOR', 'Status', [13 8]);
        setfamilydata(0, 'VCOR', 'Status', [13 9]);
        setfamilydata(0, 'CycleHCOR', 'Status', [13 8]);
        setfamilydata(0, 'CycleVCOR', 'Status', [13 9]);
        setfamilydata(0, 'BPMx', 'Status', [13 8; 13 9]);
        setfamilydata(0, 'BPMz', 'Status', [13 8; 13 9]);
    case 'S11config122' % with S11 122 BPMs
        setfamilydata(1, 'S11', 'Status')
        setfamilydata(0, 'S12', 'Status')
        setfamilydata(0, 'Q11', 'Status')
        setfamilydata(0, 'Q12', 'Status')
        setfamilydata(1, 'HCOR', 'Status', [13 8]);
        setfamilydata(1, 'VCOR', 'Status', [13 9]);
        setfamilydata(1, 'CycleHCOR', 'Status', [13 8]);
        setfamilydata(1, 'CycleVCOR', 'Status', [13 9]);
        setfamilydata(1, 'BPMx', 'Status', [13 8; 13 9]);
        setfamilydata(1, 'BPMz', 'Status', [13 8; 13 9]);
    case 'normalconfig' % without S11 120 BPMs to be obsolete
        setfamilydata(0, 'S11', 'Status')
        setfamilydata(0, 'S12', 'Status')
        setfamilydata(0, 'Q11', 'Status')
        setfamilydata(0, 'Q12', 'Status')
        setfamilydata(0, 'HCOR', 'Status', [13 8]);
        setfamilydata(0, 'VCOR', 'Status', [13 9]);
        setfamilydata(0, 'CycleHCOR', 'Status', [13 8]);
        setfamilydata(0, 'CycleVCOR', 'Status', [13 9]);
        setfamilydata(0, 'BPMx', 'Status', [13 8; 13 9]);
        setfamilydata(0, 'BPMz', 'Status', [13 8; 13 9]);
    case 'nanoscopiumconfig120' % 120 BPMs to be obsolete
        setfamilydata(1, 'S11', 'Status')
        setfamilydata(1, 'S12', 'Status')
        setfamilydata(1, 'Q11', 'Status')
        setfamilydata(1, 'Q12', 'Status')
        setfamilydata(0, 'HCOR', 'Status', [13 8]);
        setfamilydata(0, 'VCOR', 'Status', [13 9]);
        setfamilydata(0, 'CycleHCOR', 'Status', [13 8]);
        setfamilydata(0, 'CycleVCOR', 'Status', [13 9]);
        setfamilydata(0, 'BPMx', 'Status', [13 8; 13 9]);
        setfamilydata(0, 'BPMz', 'Status', [13 8; 13 9]);
    case 'nanoscopiumconfig122' % 122 BPMs 
        setfamilydata(1, 'S11', 'Status')
        setfamilydata(1, 'S12', 'Status')
        setfamilydata(1, 'Q11', 'Status')
        setfamilydata(1, 'Q12', 'Status')
        setfamilydata(0, 'HCOR', 'Status', [13 8]);
        setfamilydata(0, 'VCOR', 'Status', [13 9]);
        setfamilydata(0, 'CycleHCOR', 'Status', [13 8]);
        setfamilydata(0, 'CycleVCOR', 'Status', [13 9]);
        setfamilydata(1, 'BPMx', 'Status', [13 8; 13 9]);
        setfamilydata(1, 'BPMz', 'Status', [13 8; 13 9]);
    case 'nanoscopiumconfig122C'
        setfamilydata(1, 'S11', 'Status')
        setfamilydata(1, 'S12', 'Status')
        setfamilydata(1, 'Q11', 'Status')
        setfamilydata(1, 'Q12', 'Status')
        setfamilydata(1, 'HCOR', 'Status', [13 8]);
        setfamilydata(1, 'VCOR', 'Status', [13 9]);
        setfamilydata(1, 'CycleHCOR', 'Status', [13 8]);
        setfamilydata(1, 'CycleVCOR', 'Status', [13 9]);
        setfamilydata(1, 'BPMx', 'Status', [13 8; 13 9]);
        setfamilydata(1, 'BPMz', 'Status', [13 8; 13 9]);
    otherwise
        error('Wrong mode')
end

% switch addition corrector for HU640... TO BE REMOVED LATER
switchHU640Cor('OFF');
end
