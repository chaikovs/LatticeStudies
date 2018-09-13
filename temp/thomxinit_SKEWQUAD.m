function thomxinit(OperationalMode)
%SOLEILINIT - Initializes params for SOLEIL control in MATLAB
%

% Modified by Jianfeng Zhang @ LAL 18/06/2013
%
%==========================
% Accelerator Family Field
%==========================
% FamilyName            BPMx, HCOR, etc
% CommonNames           Shortcut name for each element
% DeviceList            [Sector, Number]
% ElementList           number in list
% Position              m, if thick, it is not the magnet center
%
% MONITOR FIELD
% Mode                  online/manual/special/simulator
% TangoNames            Device Tango Names
% Units                 Physics or HW
% HW2PhysicsFcn         function handle used to convert from hardware to physics units ==> inline will not compile, see below
% HW2PhysicsParams      params used for conversion function
% Physics2HWFcn         function handle used to convert from physics to hardware units
% Physics2HWParams      params used for conversion function
% HWUnits               units for Hardware 'A';
% PhysicsUnits          units for physics 'Rad';
% Handles               monitor handle
%
% SETPOINT FIELDS
% Mode                  online/manual/special/simulator
% TangoNames            Devices tango names
% Units                 hardware or physics
% HW2PhysicsFcn         function handle used to convert from hardware to physics units
% HW2PhysicsParams      params used for conversion function
% Physics2HWFcn         function handle used to convert from physics to hardware units
% Physics2HWParams      params used for conversion function
% HWUnits               units for Hardware 'A';
% PhysicsUnits          units for physics 'Rad';
% Range                 minsetpoint, maxsetpoint;
% Tolerance             setpoint-monitor
% Handles               setpoint handle

%=============================================
% Accelerator Toolbox Simulation Fields
%=============================================
% ATType                Quad, Sext, etc
% ATIndex               index in THERING
% ATParamGroup      param group
%
%============
% Family List
%============
%    BPMx
%    BPMz
%    HCOR
%    VCOR
%    BEND
%    Q1 to Q10
%    S1 to S10
%    RF
%    TUNE
%    DCCT
%    Machine Params
%
% NOTES
%   All sextupoles have H and V corrector and skew quadrupole windings
%
%  See Also setpathsoleil, setpathmml, aoinit, setoperationalmode, updateatindex

%
%
% TODO, Deltakick for BPM orbit response matrix  Warning optics dependent cf. Low alpha lattice
%       to be put into setoperationalmode

% DO NOT REMOVE LSN
%suppress optimization message
%#ok<*ASGLU>

% CONTROL ROOM
% Check for nanoscopium
% Check for attribute names
% Check for range value of sextupole

% If controlromm user is operator and online mode
%
%
%
%  A lot of things to be checked and modified when the ThomX machine 
%  and Tango is ready..., % Needs to be carefully tested...  by Zhang @ LAL, 02/2014.
%
%
% 24/02/2014 by Jianfeng Zhang @ LAL
%    Fix the bug to define all the family members of sextupoles.
%
%


% system gives back an visible character: carriage return!
% so comparison on the number of caracters
%  This is a very important command, decide the 
% machine mode or simualtion mode!!!!!!!!!!!!
% Need to modify it in the future... by Jianfeng Zhang @ LAL, 09/04/2014
%
[statuss WHO] = system('whoami');
if strncmp(WHO, 'operateur',9),
    ControlRoomFlag = 1;
    Mode = 'Online';
else
    ControlRoomFlag = 0;
    Mode = 'Simulator';
end

%% Default operation mode (see setoperationalmode)
if nargin < 1
    OperationalMode = 1; 
end

% Define some global variables

h = waitbar(0,'thomxinit initialization, please wait');

%==============================
%% load AcceleratorData structure
%==============================

setad([]);       %clear AcceleratorData memory
AD.SubMachine = 'StorageRing';   % Will already be defined if setpathmml was used
AD.Energy        = 50e-3; % Energy in GeV needed for magnet calibration. Do not remove!

setad(AD);

%%%%%%%%%%%%%%%%%%%%
% ACCELERATOR OBJECT
%%%%%%%%%%%%%%%%%%%%

if ~isempty(getappdata(0, 'AcceleratorObjects'))
    AO = getao;
    % Check if online and AO is from Storagering
    if ControlRoomFlag && isfield(AO, 'Q1') 
        local_tango_kill_allgroup(AO); % kill all TANGO group
    end
end
setao([]); AO =[];   %clear previous AcceleratorObjects
waitbar(0.05,h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BPM
% status field designates if BPM in use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BPMx Horizontal plane
ifam = 'BPMx';
AO.(ifam).FamilyName               = ifam;
AO.(ifam).FamilyType               = 'BPM';
AO.(ifam).MemberOf                 = {'BPM'; 'HBPM'; 'PlotFamily'; 'Archivable'};
AO.(ifam).Monitor.Mode             = Mode;
%AO.(ifam).Monitor.Mode             = 'Special';
AO.(ifam).Monitor.Units            = 'Hardware';
AO.(ifam).Monitor.HWUnits          = 'mm';
AO.(ifam).Monitor.PhysicsUnits     = 'm';
AO.(ifam).Monitor.SpecialFunctionGet = 'gethbpmgroup';
AO.(ifam).Simulated.NoiseActivated = 0; %To activate Noise on BPM reading in simulation   
%AO.(ifam).Monitor.SpecialFunctionGet = 'gethbpmaverage';
%AO.(ifam).Monitor.SpecialFunctionGet = 'gethturdevnumberyturn';

% ElemList devlist tangoname status common
% need to add the BPM inside the RF cavity
varlist = {
    1 [ 1  1] 'RI-C1/DG/BPM.010'  1 'BPMx001'
    2 [ 1  2] 'RI-C1/DG/BPM.020'  1 'BPMx002'
    3 [ 1  3] 'RI-C1/DG/BPM.030'  1 'BPMx003'
    4 [ 1  4] 'RI-C1/DG/BPM.040'  1 'BPMx004'
    5 [ 1  5] 'RI-C1/DG/BPM.050'  1 'BPMx005'
    6 [ 1  6] 'RI-C1/DG/BPM.060'  1 'BPMx006'
    7 [ 1  7] 'RI-C2/DG/BPM.010'  1 'BPMx007'
    8 [ 1  8] 'RI-C2/DG/BPM.020'  1 'BPMx008'
    9 [ 1  9] 'RI-C2/DG/BPM.030'  1 'BPMx009'
    10 [ 1  10] 'RI-C2/DG/BPM.040'  1 'BPMx010'
    11 [ 1  11] 'RI-C2/DG/BPM.050'  1 'BPMx011'
    12 [ 1  12] 'RI-C2/DG/BPM.060'  1 'BPMx012'
    };

devnumber = length(varlist);

% preallocation
AO.(ifam).ElementList = zeros(devnumber,1);
AO.(ifam).Status      = zeros(devnumber,1);
AO.(ifam).Gain        = ones(devnumber,1);
AO.(ifam).Roll        = zeros(devnumber,1);
AO.(ifam).DeviceName  = cell(devnumber,1);
AO.(ifam).CommonNames = cell(devnumber,1);
AO.(ifam).Monitor.TangoNames  = cell(devnumber,1);
AO.(ifam).Monitor.HW2PhysicsParams(:,:) = 1e-3*ones(devnumber,1);
AO.(ifam).Monitor.Physics2HWParams(:,:) = 1e3*ones(devnumber,1);
AO.(ifam).Monitor.Handles(:,1)          = NaN*ones(devnumber,1);
AO.(ifam).Monitor.DataType              = 'Scalar';

for k = 1: devnumber,
    AO.(ifam).ElementList(k)  = varlist{k,1};
    AO.(ifam).DeviceList(k,:) = varlist{k,2};
    AO.(ifam).DeviceName(k)   = deblank(varlist(k,3));
    AO.(ifam).Status(k)       = varlist{k,4};
    AO.(ifam).CommonNames(k)  = deblank(varlist(k,5));
    %need to be custumized for Tango in the future
    AO.(ifam).Monitor.TangoNames{k}  = strcat(AO.(ifam).DeviceName{k}, '/XPosSA');
end

% Group
if ControlRoomFlag
        AO.(ifam).GroupId = tango_group_create2('BPMx');
        tango_group_add(AO.(ifam).GroupId,AO.(ifam).DeviceName');
else
    AO.(ifam).GroupId = NaN;
end

AO.(ifam).History = AO.(ifam).Monitor;
AO.(ifam).History = rmfield(AO.(ifam).History,'SpecialFunctionGet');

dev = AO.(ifam).DeviceName;
AO.(ifam).History.TangoNames(:,:)       = strcat(dev, '/XPosSAHistory');
AO.(ifam).History.MemberOf = {'Plotfamily'};

% 4 electrodes?
AO.(ifam).Va = AO.(ifam).History;
AO.(ifam).Va.TangoNames(:,:)       = strcat(dev, '/VaSA');
AO.(ifam).Va.MemberOf = {'Plotfamily'};

AO.(ifam).Vb = AO.(ifam).History;
AO.(ifam).Vb.TangoNames(:,:)       = strcat(dev, '/VbSA');
AO.(ifam).Vb.MemberOf = {'Plotfamily'};

AO.(ifam).Vc = AO.(ifam).History;
AO.(ifam).Vc.TangoNames(:,:)       = strcat(dev, '/VcSA');
AO.(ifam).Vc.MemberOf = {'Plotfamily'};

AO.(ifam).Vd = AO.(ifam).History;
AO.(ifam).Vd.TangoNames(:,:)       = strcat(dev, '/VdSA');
AO.(ifam).Vd.MemberOf = {'Plotfamily'};

AO.(ifam).Sum = AO.(ifam).History;
AO.(ifam).Sum.TangoNames(:,:)       = strcat(dev, '/SumSA');
AO.(ifam).Sum.MemberOf = {'Plotfamily'};

AO.(ifam).Quad = AO.(ifam).History;
AO.(ifam).Quad.TangoNames(:,:)       = strcat(dev, '/QuadSA');
AO.(ifam).Quad.MemberOf = {'Plotfamily'};

AO.(ifam).Gaidevnumberpm = AO.(ifam).History;
AO.(ifam).Gaidevnumberpm.TangoNames(:,:)       = strcat(dev, '/Gain');
AO.(ifam).Gaidevnumberpm.MemberOf = {'Plotfamily'};

AO.(ifam).Switch = AO.(ifam).History;
AO.(ifam).Switch.TangoNames(:,:)       = strcat(dev, '/Switches');
AO.(ifam).Switch.MemberOf = {'Plotfamily'};

%% BPMz Vertical plane
ifam = 'BPMz';
AO.(ifam).FamilyName               = ifam;
AO.(ifam).FamilyType               = 'BPM';
AO.(ifam).MemberOf                 = {'BPM'; 'VBPM'; 'PlotFamily'; 'Archivable'};
AO.(ifam).Monitor.Mode             = Mode;
%AO.(ifam).Monitor.Mode             = 'Special';
AO.(ifam).Monitor.Units            = 'Hardware';
AO.(ifam).Monitor.HWUnits          = 'mm';
AO.(ifam).Monitor.PhysicsUnits     = 'm';
AO.(ifam).Monitor.SpecialFunctionGet = 'getvbpmgroup';
AO.(ifam).Simulated.NoiseActivated = 0; % To activate Noise on BPM reading in simulation     
%AO.(ifam).Monitor.SpecialFunctionGet = 'getvbpmaverage';
%AO.(ifam).Monitor.SpecialFunctionGet = 'getvturdevnumberyturn';

% devliste tangoname status common
% need to add the BPM in the RF cavity
varlist = {
    1 [ 1  1] 'RI-C1/DG/BPM.010'  1 'BPMz001'
    2 [ 1  2] 'RI-C1/DG/BPM.020'  1 'BPMz002'
    3 [ 1  3] 'RI-C1/DG/BPM.030'  1 'BPMz003'
    4 [ 1  4] 'RI-C1/DG/BPM.040'  1 'BPMz004'
    5 [ 1  5] 'RI-C1/DG/BPM.050'  1 'BPMz005'
    6 [ 1  6] 'RI-C1/DG/BPM.060'  1 'BPMz006'
    7 [ 1  7] 'RI-C2/DG/BPM.010'  1 'BPMz007'
    8 [ 1  8] 'RI-C2/DG/BPM.020'  1 'BPMz008'
    9 [ 1  9] 'RI-C2/DG/BPM.030'  1 'BPMz009'
    10 [ 1  10] 'RI-C2/DG/BPM.040'  1 'BPMz010'
    11 [ 1  11] 'RI-C2/DG/BPM.050'  1 'BPMz011'
    12 [ 1  12] 'RI-C2/DG/BPM.060'  1 'BPMz012'
    };

devnumber = length(varlist);

% preallocation
AO.(ifam).ElementList = zeros(devnumber,1);
AO.(ifam).Status      = zeros(devnumber,1);
AO.(ifam).Gain        = ones(devnumber,1);
AO.(ifam).Roll        = zeros(devnumber,1);
AO.(ifam).DeviceName  = cell(devnumber,1);
AO.(ifam).CommonNames = cell(devnumber,1);
AO.(ifam).Monitor.TangoNames  = cell(devnumber,1);
AO.(ifam).Monitor.HW2PhysicsParams(:,:) = 1e-3*ones(devnumber,1);
AO.(ifam).Monitor.Physics2HWParams(:,:) = 1e3*ones(devnumber,1);
AO.(ifam).Monitor.Handles(:,1)          = NaN*ones(devnumber,1);
AO.(ifam).Monitor.DataType              = 'Scalar';

for k = 1: devnumber,
    AO.(ifam).ElementList(k)  = varlist{k,1};
    AO.(ifam).DeviceList(k,:) = varlist{k,2};
    AO.(ifam).DeviceName(k)   = deblank(varlist(k,3));
    AO.(ifam).Status(k)       = varlist{k,4};
    AO.(ifam).CommonNames(k)  = deblank(varlist(k,5));
    AO.(ifam).Monitor.TangoNames{k}  = strcat(AO.(ifam).DeviceName{k}, '/ZPosSA');
end

% Group
if ControlRoomFlag
        AO.(ifam).GroupId = tango_group_create2('BPMz');
        tango_group_add(AO.(ifam).GroupId,AO.(ifam).DeviceName');
else
    AO.(ifam).GroupId = NaN;
end

AO.(ifam).History = AO.(ifam).Monitor;
AO.(ifam).History = rmfield(AO.(ifam).History,'SpecialFunctionGet');

dev = AO.(ifam).DeviceName;
AO.(ifam).History.TangoNames(:,:)       = strcat(dev, '/ZPosSAHistory');
AO.(ifam).History.MemberOf = {'Plotfamily'};

AO.(ifam).Va = AO.(ifam).History;
AO.(ifam).Va.TangoNames(:,:)       = strcat(dev, '/VaSA');
AO.(ifam).Va.MemberOf = {'Plotfamily'};

AO.(ifam).Vb = AO.(ifam).History;
AO.(ifam).Vb.TangoNames(:,:)       = strcat(dev, '/VbSA');
AO.(ifam).Vb.MemberOf = {'Plotfamily'};

AO.(ifam).Vc = AO.(ifam).History;
AO.(ifam).Vc.TangoNames(:,:)       = strcat(dev, '/VcSA');
AO.(ifam).Vc.MemberOf = {'Plotfamily'};

AO.(ifam).Vd = AO.(ifam).History;
AO.(ifam).Vd.TangoNames(:,:)       = strcat(dev, '/VdSA');
AO.(ifam).Vd.MemberOf = {'Plotfamily'};

AO.(ifam).Sum = AO.(ifam).History;
AO.(ifam).Sum.TangoNames(:,:)       = strcat(dev, '/SumSA');
AO.(ifam).Sum.MemberOf = {'Plotfamily'};

AO.(ifam).Quad = AO.(ifam).History;
AO.(ifam).Quad.TangoNames(:,:)       = strcat(dev, '/QuadSA');
AO.(ifam).Quad.MemberOf = {'Plotfamily'};

AO.(ifam).Gaidevnumberpm = AO.(ifam).History;
AO.(ifam).Gaidevnumberpm.TangoNames(:,:)       = strcat(dev, '/Gain');
AO.(ifam).Gaidevnumberpm.MemberOf = {'Plotfamily'};

AO.(ifam).Switch = AO.(ifam).History;
AO.(ifam).Switch.TangoNames(:,:)       = strcat(dev, '/Switches');
AO.(ifam).Switch.MemberOf = {'Plotfamily'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SLOW HORIZONTAL CORRECTORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ifam = 'HCOR';
AO.(ifam).FamilyName               = ifam;
AO.(ifam).FamilyType               = 'COR';
AO.(ifam).MemberOf                 = {'MachineConfig'; 'HCOR'; 'COR'; 'Magnet'; 'PlotFamily'; 'Archivable'};

AO.(ifam).Monitor.Mode             = Mode;
AO.(ifam).Monitor.DataType         = 'Scalar';
AO.(ifam).Monitor.Units            = 'Hardware';
AO.(ifam).Monitor.HWUnits          = 'A';
AO.(ifam).Monitor.PhysicsUnits     = 'rad';
AO.(ifam).Monitor.HW2PhysicsFcn    = @amp2k;
AO.(ifam).Monitor.Physics2HWFcn    = @k2amp;

AO.(ifam).Setpoint.Mode             = Mode;
AO.(ifam).Setpoint.DataType         = 'Scalar';
AO.(ifam).Setpoint.Units            = 'Hardware';
AO.(ifam).Setpoint.HWUnits          = 'A';
AO.(ifam).Setpoint.PhysicsUnits     = 'rad';
AO.(ifam).Setpoint.HW2PhysicsFcn    = @amp2k;
AO.(ifam).Setpoint.Physics2HWFcn    = @k2amp;

AO.(ifam).Voltage.Mode             = Mode;
AO.(ifam).Voltage.DataType         = 'Scalar';
AO.(ifam).Voltage.Units            = 'Hardware';
AO.(ifam).Voltage.HWUnits          = 'V';
AO.(ifam).Voltage.PhysicsUnits     = 'rad';
AO.(ifam).Voltage.HW2PhysicsFcn    = @amp2k;
AO.(ifam).Voltage.Physics2HWFcn    = @k2amp;

% % Get mapping from TANGO static database
AO.(ifam).Voltage.PhysicsUnits     = 'rad';
AO.(ifam).Voltage.HW2PhysicsFcn    = @amp2k;
AO.(ifam).Voltage.Physics2HWFcn    = @k2amp;

% elemlist devlist tangoname       status  common  attR
% attW      range
% the current range needed to be update for ThomX in the future
varlist = {
     1  [ 1 1] 'RI-C1/ME/STR.010      ' 1 'HCOR001' 'current   ' 'currentPM ' [-11.0  11.0]
     2  [ 1 2] 'RI-C1/ME/STR.020      ' 1 'HCOR002' 'current   ' 'currentPM ' [-11.0  11.0]
     3  [ 1 3] 'RI-C1/ME/STR.030      ' 1 'HCOR003' 'current   ' 'currentPM ' [-11.0  11.0]
     4  [ 1 4] 'RI-C1/ME/STR.040      ' 1 'HCOR004' 'current   ' 'currentPM ' [-11.0  11.0]
     5  [ 1 5] 'RI-C1/ME/STR.050      ' 1 'HCOR005' 'current   ' 'currentPM ' [-11.0  11.0]
     6  [ 1 6] 'RI-C1/ME/STR.060      ' 1 'HCOR006' 'current   ' 'currentPM ' [-11.0  11.0]
     7  [ 2 1] 'RI-C2/ME/STR.010      ' 1 'HCOR007' 'current   ' 'currentPM ' [-11.0  11.0]
     8  [ 2 2] 'RI-C2/ME/STR.020      ' 1 'HCOR008' 'current   ' 'currentPM ' [-11.0  11.0]
     9  [ 2 3] 'RI-C2/ME/STR.030      ' 1 'HCOR009' 'current   ' 'currentPM ' [-11.0  11.0]
    10  [ 2 4] 'RI-C2/ME/STR.040      ' 1 'HCOR010' 'current   ' 'currentPM ' [-11.0  11.0]
    11  [ 2 5] 'RI-C2/ME/STR.050      ' 1 'HCOR011' 'current   ' 'currentPM ' [-11.0  11.0]
    12  [ 2 6] 'RI-C2/ME/STR.060      ' 1 'HCOR012' 'current   ' 'currentPM ' [-11.0  11.0]
          };

devnumber = length(varlist);
% preallocation; meomery allocation
AO.(ifam).ElementList = zeros(devnumber,1);
AO.(ifam).Status      = zeros(devnumber,1);
AO.(ifam).Gain        = ones(devnumber,1);
AO.(ifam).Roll        = zeros(devnumber,1);
AO.(ifam).DeviceName  = cell(devnumber,1);
AO.(ifam).DeviceName  = cell(devnumber,1);
AO.(ifam).CommonNames = cell(devnumber,1);
AO.(ifam).Monitor.TangoNames  = cell(devnumber,1);
AO.(ifam).Setpoint.TangoNames = cell(devnumber,1);
AO.(ifam).Setpoint.Range = zeros(devnumber,2);

for k = 1: devnumber,
    AO.(ifam).ElementList(k)  = varlist{k,1};
    AO.(ifam).DeviceList(k,:) = varlist{k,2};
    AO.(ifam).DeviceName(k)   = deblank(varlist(k,3));
    AO.(ifam).Status(k)       = varlist{k,4};
    AO.(ifam).CommonNames(k)  = deblank(varlist(k,5));
    AO.(ifam).Monitor.TangoNames(k)  = strcat(AO.(ifam).DeviceName{k}, '/', deblank(varlist(k,6)));
    AO.(ifam).Setpoint.TangoNames(k) = strcat(AO.(ifam).DeviceName{k}, '/', deblank(varlist(k,7)));
    AO.(ifam).Setpoint.Range(k,:)      = varlist{k,8};
    % information for getrunflag
    AO.(ifam).Setpoint.RunFlagFcn = @corgetrunflag;
    AO.(ifam).Setpoint.RampRate = 1;
end

%Load fields from datablock
% AT use the "A-coefficients" for correctors plus an offset
setao(AO);
[C, Leff, MagnetType, coefficients] = magnetcoefficients(AO.(ifam).FamilyName);
for ii = 1:devnumber,
    AO.(ifam).Monitor.HW2PhysicsParams{1}(ii,:)   = coefficients/Leff;
    AO.(ifam).Monitor.Physics2HWParams{1}(ii,:)   = coefficients/Leff;
    AO.(ifam).Setpoint.HW2PhysicsParams{1}(ii,:)  = coefficients/Leff;
    AO.(ifam).Setpoint.Physics2HWParams{1}(ii,:)  = coefficients/Leff;
end

% need to customized for ThomX?
AO.(ifam).Setpoint.Tolerance(:,:)    = 1e-2*ones(devnumber,1);
% Warning optics dependent cf. Low alpha lattice
%AO.(ifam).Setpoint.DeltaRespMat(:,:) = ones(devnumber,1)*5e-5*2; % ??? urad (half used for kicking) 
 AO.(ifam).Setpoint.DeltaRespMat(:,:) = ones(devnumber,1)*5e-5*2; % 2*5 urad (half used for kicking) % LAST *5e-6*2;
%AO.(ifam).Setpoint.DeltaRespMat(:,:) = ones(devnumber,1)*0.5e-4*1; % 2*25 urad (half used for kicking)


dummy = strcat(AO.(ifam).DeviceName,'/voltage');
%AO.(ifam).Voltage.TangoNames(:,:)     = {dummy{:},['ANS-C05/EI/L-' 'HU640/currentCHE'], 'ANS-C05/EI/L-HU640/currentCHS'}';
AO.(ifam).Voltage.TangoNames(:,:)     = {dummy{:}}';

AO.(ifam).Monitor.MemberOf      = {'PlotFamily'};
AO.(ifam).Setpoint.MemberOf     = {'PlotFamily'};
AO.(ifam).Voltage.MemberOf      = {'PlotFamily'};

% Profibus configuration sync/unsync mecanism
%AO.(ifam).Profibus.BoardNumber = int32(0);
%AO.(ifam).Profibus.Group       = int32(1);
%AO.(ifam).Profibus.DeviceName  = 'ANS/AE/DP.CH';

% Group
if ControlRoomFlag
        AO.(ifam).GroupId = tango_group_create2(ifam);
        tango_group_add(AO.(ifam).GroupId,AO.(ifam).DeviceName(find(AO.(ifam).Status),:)');
else
    AO.(ifam).GroupId = nan;
end

%convert response matrix kicks to HWUnits (after AO is loaded to AppData)
setao(AO);   %required to make physics2hw function
AO.(ifam).Setpoint.DeltaRespMat = physics2hw(AO.(ifam).FamilyName,'Setpoint', ...
    AO.(ifam).Setpoint.DeltaRespMat, AO.(ifam).DeviceList);

AO.(ifam).TangoSetpoint = AO.(ifam).Setpoint;
AO.(ifam).TangoSetpoint.TangoNames(:,:)    = ...
    strcat(AO.(ifam).DeviceName,'/currentTotal');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SLOW VERTICAL CORRECTORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ifam = 'VCOR';
AO.(ifam).FamilyName               = ifam;
AO.(ifam).FamilyType               = 'COR';
AO.(ifam).MemberOf                 = {'MachineConfig'; 'VCOR'; 'COR'; 'Magnet'; 'PlotFamily'; 'Archivable'};

AO.(ifam).Monitor.Mode             = Mode;
AO.(ifam).Monitor.DataType         = 'Scalar';
AO.(ifam).Monitor.Units            = 'Hardware';
AO.(ifam).Monitor.HWUnits          = 'A';
AO.(ifam).Monitor.PhysicsUnits     = 'meter^-2';
AO.(ifam).Monitor.HW2PhysicsFcn = @amp2k;
AO.(ifam).Monitor.Physics2HWFcn = @k2amp;

AO.(ifam).Setpoint.Mode             = Mode;
AO.(ifam).Setpoint.DataType         = 'Scalar';
AO.(ifam).Setpoint.Units            = 'Hardware';
AO.(ifam).Setpoint.HWUnits          = 'A';
AO.(ifam).Setpoint.PhysicsUnits     = 'rad';
AO.(ifam).Setpoint.HW2PhysicsFcn    = @amp2k;
AO.(ifam).Setpoint.Physics2HWFcn    = @k2amp;

AO.(ifam).Voltage.Mode             = Mode;
AO.(ifam).Voltage.DataType         = 'Scalar';
AO.(ifam).Voltage.Units            = 'Hardware';
AO.(ifam).Voltage.HWUnits          = 'A';
AO.(ifam).Voltage.PhysicsUnits     = 'rad';
AO.(ifam).Voltage.HW2PhysicsFcn    = @amp2k;
AO.(ifam).Voltage.Physics2HWFcn    = @k2amp;

varlist = {
     1  [ 1 1] 'RI-C1/ME/STR.010      ' 1 'VCOR001' 'current   ' 'currentPM ' [-11.0  11.0]
     2  [ 1 2] 'RI-C1/ME/STR.020      ' 1 'VCOR002' 'current   ' 'currentPM ' [-11.0  11.0]
     3  [ 1 3] 'RI-C1/ME/STR.030      ' 1 'VCOR003' 'current   ' 'currentPM ' [-11.0  11.0]
     4  [ 1 4] 'RI-C1/ME/STR.040      ' 1 'VCOR004' 'current   ' 'currentPM ' [-11.0  11.0]
     5  [ 1 5] 'RI-C1/ME/STR.050      ' 1 'VCOR005' 'current   ' 'currentPM ' [-11.0  11.0]
     6  [ 1 6] 'RI-C1/ME/STR.060      ' 1 'VCOR006' 'current   ' 'currentPM ' [-11.0  11.0]
     7  [ 2 1] 'RI-C2/ME/STR.010      ' 1 'VCOR007' 'current   ' 'currentPM ' [-11.0  11.0]
     8  [ 2 2] 'RI-C2/ME/STR.020      ' 1 'VCOR008' 'current   ' 'currentPM ' [-11.0  11.0]
     9  [ 2 3] 'RI-C2/ME/STR.030      ' 1 'VCOR009' 'current   ' 'currentPM ' [-11.0  11.0]
    10  [ 2 4] 'RI-C2/ME/STR.040      ' 1 'VCOR010' 'current   ' 'currentPM ' [-11.0  11.0]
    11  [ 2 5] 'RI-C2/ME/STR.050      ' 1 'VCOR011' 'current   ' 'currentPM ' [-11.0  11.0]
    12  [ 2 6] 'RI-C2/ME/STR.060      ' 1 'VCOR012' 'current   ' 'currentPM ' [-11.0  11.0]
          };

devnumber = length(varlist);
% preallocation
AO.(ifam).ElementList = zeros(devnumber,1);
AO.(ifam).Status      = zeros(devnumber,1);
AO.(ifam).Gain        = ones(devnumber,1);
AO.(ifam).Roll        = zeros(devnumber,1);
AO.(ifam).DeviceName  = cell(devnumber,1);
AO.(ifam).DeviceName  = cell(devnumber,1);
AO.(ifam).CommonNames = cell(devnumber,1);
AO.(ifam).Monitor.TangoNames  = cell(devnumber,1);
AO.(ifam).Setpoint.TangoNames = cell(devnumber,1);
AO.(ifam).Setpoint.Range = zeros(devnumber,2);

for k = 1: devnumber,
    AO.(ifam).ElementList(k)  = varlist{k,1};
    AO.(ifam).DeviceList(k,:) = varlist{k,2};
    AO.(ifam).DeviceName(k)   = deblank(varlist(k,3));
    AO.(ifam).Status(k)       = varlist{k,4};
    AO.(ifam).CommonNames(k)  = deblank(varlist(k,5));
    AO.(ifam).Monitor.TangoNames(k)  = strcat(AO.(ifam).DeviceName{k}, '/', deblank(varlist(k,6)));
    AO.(ifam).Setpoint.TangoNames(k) = strcat(AO.(ifam).DeviceName{k}, '/', deblank(varlist(k,7)));
    AO.(ifam).Setpoint.Range(k,:)      = varlist{k,8};
    
    % information for getrunflag
    AO.(ifam).Setpoint.RunFlagFcn = @corgetrunflag;
    AO.(ifam).Setpoint.RampRate = 1;
end

%Load fields from datablock
% AT use the "A-coefficients" for correctors plus an offset
setao(AO);
[C, Leff, MagnetType, coefficients] = magnetcoefficients(AO.(ifam).FamilyName);

for ii = 1:devnumber,
    AO.(ifam).Monitor.HW2PhysicsParams{1}(ii,:)   = coefficients/Leff;
    AO.(ifam).Monitor.Physics2HWParams{1}(ii,:)   = coefficients/Leff;
    AO.(ifam).Setpoint.HW2PhysicsParams{1}(ii,:)  = coefficients/Leff;
    AO.(ifam).Setpoint.Physics2HWParams{1}(ii,:)  = coefficients/Leff;
end

AO.(ifam).Monitor.MemberOf  = {'PlotFamily'};
AO.(ifam).Voltage.MemberOf  = {'PlotFamily'};
AO.(ifam).Setpoint.MemberOf  = {'PlotFamily'};

AO.(ifam).Setpoint.Tolerance(:,:)    = 1e-2*ones(devnumber,1);
% Warning optics dependent cf. Low alpha lattice
%AO.(ifam).Setpoint.DeltaRespMat(:,:) = ones(devnumber,1)*10e-5*2;
 AO.(ifam).Setpoint.DeltaRespMat(:,:) = ones(devnumber,1)*5e-5*2; % LAST *10e-6*2;

dummy = strcat(AO.(ifam).DeviceName,'/voltage');
AO.(ifam).Voltage.TangoNames(:,:)     = {dummy{:}};

% Profibus configuration sync/yn=unsync mecanism
%AO.(ifam).Profibus.BoardNumber = int32(0);
%AO.(ifam).Profibus.Group       = int32(1);
%AO.(ifam).Profibus.DeviceName  = 'ANS/AE/DP.CV';


% Group
if ControlRoomFlag
        AO.(ifam).GroupId = tango_group_create2(ifam);
        tango_group_add(AO.(ifam).GroupId,AO.(ifam).DeviceName(find(AO.(ifam).Status),:)');
else
    AO.(ifam).GroupId = nan;
end

%convert response matrix kicks to HWUnits (after AO is loaded to AppData)
setao(AO);   %required to make physics2hw function

AO.(ifam).Setpoint.DeltaRespMat = physics2hw(AO.(ifam).FamilyName,'Setpoint', ...
    AO.(ifam).Setpoint.DeltaRespMat, AO.(ifam).DeviceList);
AO.(ifam).TangoSetpoint = AO.(ifam).Setpoint;
AO.(ifam).TangoSetpoint.TangoNames(:,:)    = ...
    strcat(AO.(ifam).DeviceName,'/currentTotal');


%=============================
%        MAIN MAGNETS
%
%   Need to customize for ThomX current range
%=============================

%==============
%% DIPOLES
%==============

varlist={
     1  [ 1  1] 'RI-C1/ME/DP.010-DC' 1 'BEND.01' [+0 +560]
     2  [ 1  2] 'RI-C1/ME/DP.020-DC' 1 'BEND.02' [+0 +560]
     3  [ 1  3] 'RI-C1/ME/DP.030-DC' 1 'BEND.03' [+0 +560]
     4  [ 1  4] 'RI-C1/ME/DP.040-DC' 1 'BEND.04' [+0 +560]
     5  [ 2  1] 'RI-C2/ME/DP.010-DC' 1 'BEND.05' [+0 +560]
     6  [ 2  2] 'RI-C2/ME/DP.020-DC' 1 'BEND.06' [+0 +560]
     7  [ 2  3] 'RI-C2/ME/DP.030-DC' 1 'BEND.07' [+0 +560]
     8  [ 2  4] 'RI-C2/ME/DP.040-DC' 1 'BEND.08' [+0 +560]
};

% *** BEND ***
ifam='BEND';
AO.(ifam).FamilyName                 = ifam;
AO.(ifam).FamilyType                 = 'BEND';
AO.(ifam).MemberOf                   = {'MachineConfig'; 'BEND'; 'Magnet'; 'PlotFamily'; 'Archivable'};
HW2PhysicsParams                    = magnetcoefficients('BEND');
Physics2HWParams                    = HW2PhysicsParams;

AO.(ifam).Monitor.Mode               = Mode;
AO.(ifam).Monitor.DataType           = 'Scalar';
AO.(ifam).Monitor.Units              = 'Hardware';
AO.(ifam).Monitor.HW2PhysicsFcn      = @bend2gev;
AO.(ifam).Monitor.Physics2HWFcn      = @gev2bend;
AO.(ifam).Monitor.HWUnits            = 'A';
AO.(ifam).Monitor.PhysicsUnits       = 'GeV';

devnumber = size(varlist,1);
% preallocation
AO.(ifam).ElementList         = zeros(devnumber,1);
AO.(ifam).Status              = zeros(devnumber,1);
AO.(ifam).DeviceName          = cell(devnumber,1);
AO.(ifam).DeviceName          = cell(devnumber,1);
AO.(ifam).CommonNames         = cell(devnumber,1);
AO.(ifam).Monitor.TangoNames  = cell(devnumber,1);
% make Setpoint structure than specific data
AO.(ifam).Setpoint            = AO.(ifam).Monitor;
AO.(ifam).Setpoint.Range      = zeros(devnumber,2);

for ik = 1: devnumber,   
    AO.(ifam).ElementList(ik)           = varlist{ik,1};
    AO.(ifam).DeviceList(ik,:)          = varlist{ik,2};
    AO.(ifam).DeviceName(ik)            = deblank(varlist(ik,3));
    AO.(ifam).Status(ik)                = varlist{ik,4};
    AO.(ifam).CommonNames(ik)           = deblank(varlist(ik,5));
    AO.(ifam).Monitor.TangoNames(ik)    = strcat(AO.(ifam).DeviceName{ik}, {'/current'});   
    AO.(ifam).Setpoint.TangoNames(ik)   = strcat(AO.(ifam).DeviceName{ik}, {'/currentPM'});
    AO.(ifam).Status(ik)                = varlist{ik,4};
    AO.(ifam).Setpoint.Range(ik,:)      = varlist{ik,6};
end

AO.(ifam).Monitor.Handles(:,1) = NaN*ones(devnumber,1);

HW2PhysicsParams = magnetcoefficients(AO.(ifam).FamilyName);
Physics2HWParams = magnetcoefficients(AO.(ifam).FamilyName);

val = 1;
for ii=1:devnumber,
    AO.(ifam).Monitor.HW2PhysicsParams{1}(ii,:)                 = HW2PhysicsParams;
    AO.(ifam).Monitor.HW2PhysicsParams{2}(ii,:)                 = val;
    AO.(ifam).Monitor.Physics2HWParams{1}(ii,:)                 = Physics2HWParams;
    AO.(ifam).Monitor.Physics2HWParams{2}(ii,:)                 = val;
end
% same configuration for Monitor and Setpoint value concerning hardware to physics units
AO.(ifam).Setpoint.HW2PhysicsParams = AO.(ifam).Monitor.HW2PhysicsParams;
AO.(ifam).Setpoint.Physics2HWParams = AO.(ifam).Monitor.Physics2HWParams;

AO.(ifam).Setpoint = AO.(ifam).Monitor;
AO.(ifam).Desired  = AO.(ifam).Monitor;
AO.(ifam).Setpoint.MemberOf  = {'PlotFamily'};
AO.(ifam).Setpoint.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/currentPM');

AO.(ifam).Setpoint.Tolerance(:,:) = 0.05;
AO.(ifam).Setpoint.DeltaRespMat(:,:) = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QUADRUPOLE MAGNETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear varlist;
varlist.QP1={
    1 [ 1  1 ] 'RI-C1/ME/QP.010   ' 1 ' QP1.1' [-250 +0]
    2 [ 1  12] 'RI-C1/ME/QP.120   ' 1 ' QP1.2' [-250 +0]
    3 [ 2  1 ] 'RI-C2/ME/QP.010   ' 1 ' QP1.3' [-250 +0]
    4 [ 2  12] 'RI-C2/ME/QP.120   ' 1 ' QP1.4' [-250 +0]
    };

varlist.QP2={
    1 [ 1  2 ] 'RI-C1/ME/QP.020   ' 1 ' QP2.1' [+0 +250]
    2 [ 1  11] 'RI-C1/ME/QP.110   ' 1 ' QP2.2' [+0 +250]
    3 [ 2  2 ] 'RI-C2/ME/QP.020   ' 1 ' QP2.3' [+0 +250]
    4 [ 2  11] 'RI-C2/ME/QP.110   ' 1 ' QP2.4' [+0 +250]
    };

varlist.QP31={
    1 [ 1  3 ] 'RI-C1/ME/QP.030   ' 1 ' QP31.1' [-250 +0]
    2 [ 1  10] 'RI-C1/ME/QP.100   ' 1 ' QP31.2' [-250 +0]
    3 [ 2  3 ] 'RI-C2/ME/QP.030   ' 1 ' QP31.3' [-250 +0]
    4 [ 2  10] 'RI-C2/ME/QP.100   ' 1 ' QP31.4' [-250 +0]
    };

varlist.QP41={
    1 [ 1  4 ] 'RI-C1/ME/QP.040   ' 1 ' QP41.1' [+0 +250]
    2 [ 1  9 ] 'RI-C1/ME/QP.090   ' 1 ' QP41.2' [+0 +250]
    3 [ 2  4 ] 'RI-C2/ME/QP.040   ' 1 ' QP41.3' [+0 +250]
    4 [ 2  9 ] 'RI-C2/ME/QP.090   ' 1 ' QP41.4' [+0 +250]
    };

varlist.QP4={
    1 [ 1  5 ] 'RI-C1/ME/QP.050   ' 1 ' QP4.1' [+0 +250]
    2 [ 1  8 ] 'RI-C1/ME/QP.080   ' 1 ' QP4.2' [+0 +250]
    3 [ 2  5 ] 'RI-C2/ME/QP.050   ' 1 ' QP4.3' [+0 +250]
    4 [ 2  8 ] 'RI-C2/ME/QP.080   ' 1 ' QP4.4' [+0 +250]
    };

varlist.QP3={
    1 [ 1  6 ] 'RI-C1/ME/QP.060   ' 1 ' QP3.1' [-250 +0]
    2 [ 1  7 ] 'RI-C1/ME/QP.070   ' 1 ' QP3.2' [-250 +0]
    3 [ 2  6 ] 'RI-C2/ME/QP.060   ' 1 ' QP3.3' [-250 +0]
    4 [ 2  7 ] 'RI-C2/ME/QP.070   ' 1 ' QP3.4' [-250 +0]
    };





for k = 1:6,

  % set up for QP1, QP2, QP3, QP4
  ifam = ['QP' num2str(k)];
  %set up for QP31, QP41
  if (k==5)    ifam = 'QP31'; end;
  if (k==6)    ifam = 'QP41'; end;
    
    AO.(ifam).FamilyName                 = ifam;
    AO.(ifam).FamilyType                 = 'QUAD';
    AO.(ifam).MemberOf                   = {'MachineConfig'; 'QUAD'; 'Magnet'; 'PlotFamily'; 'Archivable'};
    AO.(ifam).Monitor.Mode               = Mode;
    AO.(ifam).Monitor.DataType           = 'Scalar';
    AO.(ifam).Monitor.Units              = 'Hardware';
    AO.(ifam).Monitor.HWUnits            = 'A';
    AO.(ifam).Monitor.PhysicsUnits       = 'meter^-2';
    AO.(ifam).Monitor.HW2PhysicsFcn      = @amp2k;
    AO.(ifam).Monitor.Physics2HWFcn      = @k2amp;
    
    %ifam = sprintf('QP%s', num2str(k));

    % set the parameters of each device in the family
    devnumber = size(varlist.(ifam),1);
    % preallocation
    AO.(ifam).ElementList = zeros(devnumber,1);
    AO.(ifam).Status      = zeros(devnumber,1);
    AO.(ifam).DeviceName  = cell(devnumber,1);
    AO.(ifam).CommonNames = cell(devnumber,1);
    AO.(ifam).Monitor.TangoNames  = cell(devnumber,1);
    % make Setpoint structure than specific data
    AO.(ifam).Setpoint = AO.(ifam).Monitor;
    AO.(ifam).Setpoint.Range = zeros(devnumber,2);
    for ik = 1: devnumber,
        AO.(ifam).ElementList(ik)  = varlist.(ifam){ik,1};
        AO.(ifam).DeviceList(ik,:) = varlist.(ifam){ik,2};
        AO.(ifam).DeviceName(ik)   = deblank(varlist.(ifam)(ik,3));
        AO.(ifam).Status(ik)       = varlist.(ifam){ik,4};
        AO.(ifam).CommonNames(ik)  = deblank(varlist.(ifam)(ik,5));
        AO.(ifam).Monitor.TangoNames(ik)  = strcat(AO.(ifam).DeviceName{ik}, {'/currentPM'});
        AO.(ifam).Setpoint.TangoNames(ik) = strcat(AO.(ifam).DeviceName{ik}, {'/currentPM'});
        AO.(ifam).Setpoint.Range(ik,:)      = varlist.(ifam){ik,6};
    end
    
    %AO.(ifam).Monitor.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,{'/currentPM'});
    AO.(ifam).Monitor.Handles(:,1) = NaN*ones(devnumber,1);
   
    % Conversion between HW and Physics
    HW2PhysicsParams = magnetcoefficients(AO.(ifam).FamilyName);
    Physics2HWParams = magnetcoefficients(AO.(ifam).FamilyName);
    
    val = 1.0; % scaling factor
    for ii=1:devnumber,
        AO.(ifam).Monitor.HW2PhysicsParams{1}(ii,:)                 = HW2PhysicsParams';
        AO.(ifam).Monitor.HW2PhysicsParams{2}(ii,:)                 = val;
        AO.(ifam).Monitor.Physics2HWParams{1}(ii,:)                 = Physics2HWParams';
        AO.(ifam).Monitor.Physics2HWParams{2}(ii,:)                 = val;
    end
    % same configuration for Monitor and Setpoint value concerning hardware to physics units
    AO.(ifam).Setpoint.HW2PhysicsParams = AO.(ifam).Monitor.HW2PhysicsParams;
    AO.(ifam).Setpoint.Physics2HWParams = AO.(ifam).Monitor.Physics2HWParams;
    
    % Group
    if ControlRoomFlag
            AO.(ifam).GroupId = tango_group_create2(ifam);
            tango_group_add(AO.(ifam).GroupId,AO.(ifam).DeviceName');
    else
        AO.(ifam).GroupId = nan;
    end
    
    % to be part of plotfamily
    AO.(ifam).Setpoint.MemberOf  = {'PlotFamily'};
    % set tolerance for setting values
    AO.(ifam).Setpoint.Tolerance(:,:) = 0.02*ones(devnumber,1);
    % information for getrunflag
    AO.(ifam).Setpoint.RunFlagFcn = @tangogetrunflag;
    AO.(ifam).Setpoint.RampRate = 1;
       
    AO.(ifam).Desired  = AO.(ifam).Monitor;
    
    AO.(ifam).Voltage = AO.(ifam).Monitor;
    AO.(ifam).Voltage.MemberOf  = {'PlotFamily'};
    AO.(ifam).Voltage.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/voltage');
    AO.(ifam).Voltage.HWUnits            = 'V';
    
    AO.(ifam).DCCT = AO.(ifam).Monitor;
    AO.(ifam).DCCT.MemberOf  = {'PlotFamily'};
    AO.(ifam).DCCT.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/current');
    
    AO.(ifam).CurrentTotal = AO.(ifam).Monitor;
    AO.(ifam).CurrentTotal.MemberOf  = {'PlotFamily'};
    AO.(ifam).CurrentTotal.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/currentTotal');
    
    AO.(ifam).CurrentSetpoint = AO.(ifam).Monitor;
    AO.(ifam).CurrentSetpoint.MemberOf  = {'PlotFamily'};
    AO.(ifam).CurrentSetpoint.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/currentSetpoint');
    
    AO.(ifam).SumOffset = AO.(ifam).Monitor;
    AO.(ifam).SumOffset.MemberOf  = {'PlotFamily'};
    AO.(ifam).SumOffset.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/ecart3');
    
    AO.(ifam).Offset1 = AO.(ifam).Monitor;
    AO.(ifam).Offset1.MemberOf  = {'PlotFamily'};
    AO.(ifam).Offset1.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/currentOffset1');
     
    
    %  % Profibus configuration
% need to customized for ThomX
    % if k < 6
    %     AO.(ifam).Profibus.BoardNumber = int32(1);
    %     AO.(ifam).Profibus.Group       = int32(k);
    %     AO.(ifam).Profibus.DeviceName  = 'ANS/AE/DP.QP';
    % elseif k < 11
    %     AO.(ifam).Profibus.BoardNumber = int32(0);
    %     AO.(ifam).Profibus.Group       = int32(k-5);
    %     AO.(ifam).Profibus.DeviceName  = 'ANS/AE/DP.QP';
    % else % NANOSCOPIUM
    %     AO.(ifam).Profibus.BoardNumber = int32(0);
    %     AO.(ifam).Profibus.Group       = int32(6);
    %     AO.(ifam).Profibus.DeviceName  = 'ANS-C13/AE/DP.NANOSCOPIUM';
    % end

    
        % TODO : this is lattice dependent !
% step size of the quadrupole strength to get the tune response
% matrix do tune correction, etc. (Used by steptune.m)
% The current value is SOLEIL values, need
% to customized for THOMX!!!!
% step for tuneshift of 1-e-2 in one of the planes
% Called by measrespmat.m

%now use the the physical unit , 1 A = 1k
% Need to change to hard ware value in the future
% when the ThomX magnet is ready
deltak = 1e-6; % physics unit = hardware unit, need to modify for the real machine in the future
AO.(ifam).Setpoint.DeltaRespMat = deltak; % T/m^-2 = A


    % deltaI [ampere] change of the quadrupole for betatron
    % function measurement (Used by measbeta.m).
    % need to customize with the HW units for ThomX in the future
    % when the machine is ready
    
    AO.(ifam).Setpoint.DeltaKBeta = 0.005;  % maximum physics value, otherwise, the 
                                           % linear equation is not valid
                                           % anymore!!!!!


            %convert response matrix kicks to HWUnits (after AO is loaded to AppData)
    setao(AO);   %required to make physics2hw function

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Skew Quadrupole data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ifam = 'QT';
AO.(ifam).FamilyName               = ifam;
AO.(ifam).FamilyType               = 'SkewQuad';
AO.(ifam).MemberOf                 = {'MachineConfig'; 'Magnet'; 'PlotFamily'; 'Archivable'};

AO.(ifam).Monitor.Mode             = Mode;
AO.(ifam).Monitor.DataType         = 'Scalar';
AO.(ifam).Monitor.Units            = 'Hardware';
AO.(ifam).Monitor.HWUnits          = 'A';
AO.(ifam).Monitor.PhysicsUnits     = 'meter^-2';
AO.(ifam).Monitor.HW2PhysicsFcn = @amp2k;
AO.(ifam).Monitor.Physics2HWFcn = @k2amp;

% elemlist devlist tangoname status  common   range
clear varlist;
varlist = {
     1  [ 1 1] 'RI-C1/ME/QT.010          ' 1 'QT001'  [ -7.0   7.0]
     2  [ 1 2] 'RI-C1/ME/QT.020          ' 1 'QT002'  [ -7.0   7.0]
     3  [ 1 3] 'RI-C1/ME/QT.030          ' 1 'QT003'  [ -7.0   7.0]
     4  [ 1 4] 'RI-C1/ME/QT.040          ' 1 'QT004'  [ -7.0   7.0]
     5  [2 1] 'RI-C2/ME/QT.050           ' 1 'QT005'  [ -7.0   7.0]
     6  [2 2] 'RI-C2/ME/QT.060           ' 1 'QT006'  [ -7.0   7.0]
     7  [2 3] 'RI-C2/ME/QT.070           ' 1 'QT007'  [ -7.0   7.0]
     8  [2 4] 'RI-C2/ME/QT.080           ' 1 'QT008'  [ -7.0   7.0]
    };

devnumber = size(varlist,1);
% preallocation
AO.(ifam).ElementList = zeros(devnumber,1);
AO.(ifam).Status      = zeros(devnumber,1);
AO.(ifam).Gain        = ones(devnumber,1);
AO.(ifam).Roll        = zeros(devnumber,1);
AO.(ifam).DeviceName  = cell(devnumber,1);
AO.(ifam).DeviceName  = cell(devnumber,1);
AO.(ifam).CommonNames = cell(devnumber,1);
AO.(ifam).Setpoint = AO.(ifam).Monitor;
AO.(ifam).Monitor.TangoNames  = cell(devnumber,1);
AO.(ifam).Setpoint.TangoNames = cell(devnumber,1);
AO.(ifam).Setpoint.Range = zeros(devnumber,2);
AO.(ifam).Monitor.Handles(:,1)    = NaN*ones(devnumber,1);

for k = 1: devnumber,
    AO.(ifam).ElementList(k)  = varlist{k,1};
    AO.(ifam).DeviceList(k,:) = varlist{k,2};
    AO.(ifam).DeviceName(k)   = deblank(varlist(k,3));
    AO.(ifam).Status(k)       = varlist{k,4};
    AO.(ifam).CommonNames(k)  = deblank(varlist(k,5));
    AO.(ifam).Monitor.TangoNames{k}  = strcat(AO.(ifam).DeviceName{k}, '/currentPM');
    AO.(ifam).Setpoint.TangoNames{k} = strcat(AO.(ifam).DeviceName{k}, '/currentPM');
    AO.(ifam).Setpoint.Range(k,:)      = varlist{k,6};
    % information for getrunflag
    AO.(ifam).Setpoint.RunFlagFcn = @tangogetrunflag;
    AO.(ifam).Setpoint.RampRate = 1;
end
% Load coeeficients fot thin element
coefficients = magnetcoefficients(AO.(ifam).FamilyName);

for ii=1:devnumber,
    AO.(ifam).Monitor.HW2PhysicsParams{1}(ii,:)  = coefficients;
    AO.(ifam).Monitor.Physics2HWParams{1}(ii,:)  = coefficients;
    AO.(ifam).Setpoint.HW2PhysicsParams{1}(ii,:)  = coefficients;
    AO.(ifam).Setpoint.Physics2HWParams{1}(ii,:)  = coefficients;
end

AO.(ifam).Desired = AO.(ifam).Monitor;
AO.(ifam).Setpoint.MemberOf  = {'PlotFamily'};
AO.(ifam).Setpoint.TangoNames(:,:)    = strcat(AO.(ifam).DeviceName,'/currentPM');

AO.(ifam).Setpoint.Tolerance(:,:) = 1000*ones(devnumber,1);
AO.(ifam).Setpoint.DeltaRespMat(:,:) = 3*ones(devnumber,1);
AO.(ifam).Setpoint.DeltaSkewK = 1; % for SkewQuad efficiency toward dispersion .
% information for getrunflag
AO.(ifam).Setpoint.RunFlagFcn = @tangogetrunflag;
AO.(ifam).Setpoint.RampRate = 1;

% Profibus configuration
% AO.(ifam).Profibus.BoardNumber = int32(0);
% AO.(ifam).Profibus.Group       = int32(1);
% AO.(ifam).Profibus.DeviceName  = 'ANS/AE/DP.QT';

% Group
if ControlRoomFlag
        AO.(ifam).GroupId = tango_group_create2(ifam);
        tango_group_add(AO.(ifam).GroupId,AO.(ifam).DeviceName(find(AO.(ifam).Status),:)');
else
    AO.(ifam).GroupId = nan;
end

%convert response matrix kicks to HWUnits (after AO is loaded to AppData)
setao(AO);   %required to make physics2hw function
AO.(ifam).Setpoint.DeltaRespMat = physics2hw(AO.(ifam).FamilyName, 'Setpoint', ...
    AO.(ifam).Setpoint.DeltaRespMat, AO.(ifam).DeviceList);



%% All quadrupoles
ifam = 'Qall';
AO.(ifam).FamilyName             = ifam;
AO.(ifam).MemberOf               = {ifam;'PlotFamily'};
AO.(ifam).FamilyType             = 'Qall';
AO.(ifam).Mode                   = Mode;

AO.(ifam).DeviceName= {};

%%generate the element list of all the quadrupoles
% configurations for the "plotfamily" feature
for k = 1:6, 
    iquad = sprintf('QP%d',k);
%set up for QP31, QP41
  if (k==5)    iquad = 'QP31'; end;
  if (k==6)    iquad = 'QP41'; end;
    
    AO.(ifam).DeviceName = {AO.(ifam).DeviceName{:} AO.(iquad).DeviceName{:}};
end

AO.(ifam).DeviceName = AO.(ifam).DeviceName';

AO.(ifam).DeviceList = [];

devnumber = length(AO.(ifam).DeviceName);

% build fake device list in order to use plotfamily
for k =1:devnumber,
    AO.(ifam).DeviceList = [AO.(ifam).DeviceList; 1 k];
end

AO.(ifam).Monitor.HW2PhysicsParams(:,:) = 1e-3*ones(devnumber,1);
AO.(ifam).Monitor.Physics2HWParams(:,:) = 1e3*ones(devnumber,1);

AO.(ifam).ElementList            = (1:devnumber)';
AO.(ifam).Status                 = ones(devnumber,1);
AO.(ifam).Monitor.Mode           = Mode;
AO.(ifam).Monitor.DataType       = 'Scalar';
AO.(ifam).Monitor.Units          = 'Hardware';
AO.(ifam).Monitor.HWUnits        = 'ampere';
AO.(ifam).Monitor.PhysicsUnits   = 'meter^-2';
AO.(ifam).Monitor.TangoNames     = strcat(AO.(ifam).DeviceName, '/current');
AO.(ifam).Setpoint               = AO.(ifam).Monitor;
AO.(ifam).Setpoint.MemberOf      = {'PlotFamily'};
AO.(ifam).Setpoint.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/currentPM');

AO.(ifam).Voltage = AO.(ifam).Monitor;
AO.(ifam).Voltage.MemberOf  = {'PlotFamily'};
AO.(ifam).Voltage.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/voltage');

AO.(ifam).DCCT = AO.(ifam).Monitor;
AO.(ifam).DCCT.MemberOf  = {'PlotFamily'};
AO.(ifam).DCCT.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/current');

AO.(ifam).CurrentTotal = AO.(ifam).Monitor;
AO.(ifam).CurrentTotal.MemberOf  = {'PlotFamily'};
AO.(ifam).CurrentTotal.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/currentTotal');

AO.(ifam).CurrentSetpoint = AO.(ifam).Monitor;
AO.(ifam).CurrentSetpoint.MemberOf  = {'PlotFamily'};
AO.(ifam).CurrentSetpoint.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/currentSetpoint');

AO.(ifam).SumOffset = AO.(ifam).Monitor;
AO.(ifam).SumOffset.MemberOf  = {'PlotFamily'};
AO.(ifam).SumOffset.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/ecart3');

AO.(ifam).Offset1 = AO.(ifam).Monitor;
AO.(ifam).Offset1.MemberOf  = {'PlotFamily'};
AO.(ifam).Offset1.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/currentOffset1');

clear tune
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SEXTUPOLE MAGNETS; maybe there are some problems 
%   in this section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear varlist;
varlist.SX1 ={
    1 [ 1  1 ] 'RI-C1/ME/SP.010  ' 1 'SX1.1' [-349 +000]
    2 [ 1  6 ] 'RI-C1/ME/SP.060  ' 1 'SX1.2' [-349 +000]
    3 [ 2  1 ] 'RI-C2/ME/SP.010  ' 1 'SX1.3' [-349 +000]
    4 [ 2  6 ] 'RI-C2/ME/SP.060  ' 1 'SX1.4' [-349 +000]
};


varlist.SX2 ={
    1 [ 1  3 ] 'RI-C1/ME/SP.030  ' 1 'SX2.1' [+000 +349]
    2 [ 1  4 ] 'RI-C1/ME/SP.040  ' 1 'SX2.2' [+000 +349]
    3 [ 2  3 ] 'RI-C2/ME/SP.030  ' 1 'SX2.3' [+000 +349]
    4 [ 2  4 ] 'RI-C2/ME/SP.040  ' 1 'SX2.4' [+000 +349]
};


varlist.SX3 ={
    1 [ 1  2 ] 'RI-C1/ME/SP.020  ' 1 'SX3.1' [-349 +000]
    2 [ 1  5 ] 'RI-C1/ME/SP.050  ' 1 'SX3.2' [-349 +000]
    3 [ 2  2 ] 'RI-C2/ME/SP.020  ' 1 'SX3.3' [-349 +000]
    4 [ 2  5 ] 'RI-C2/ME/SP.050  ' 1 'SX3.4' [-349 +000]
};

% varlist={
%     1 [ 1  1 ] 'RI-C1/ME/SP.010  ' 1 'SX1' [-349 +000]
%     1 [ 1  3 ] 'RI-C1/ME/SP.030  ' 1 'SX2' [+000 +349]
%     1 [ 1  2 ] 'RI-C1/ME/SP.020  ' 1 'SX3' [-349 +000]
%         };

for k = 1:3,
    ifam = ['SX' num2str(k)];
   %  ifam = cell2mat(deblank(varlist(k,5)));
     
    AO.(ifam).FamilyName                = ifam;
    AO.(ifam).FamilyType                = 'SEXT';
    AO.(ifam).MemberOf                  = {'MachineConfig'; 'SEXT'; 'Magnet';'PlotFamily'; 'Archivable'};
    AO.(ifam).Monitor.Mode              = Mode;
    AO.(ifam).Monitor.DataType          = 'Scalar';
    AO.(ifam).Monitor.Units             = 'Hardware';
    AO.(ifam).Monitor.HWUnits           = 'A';
    AO.(ifam).Monitor.PhysicsUnits      = 'meter^-3';
    AO.(ifam).Monitor.HW2PhysicsFcn     = @amp2k;
    AO.(ifam).Monitor.Physics2HWFcn     = @k2amp;
    
    
    % set the parameters of each device in the family
    devnumber = size(varlist.(ifam),1);
    % preallocation
    AO.(ifam).ElementList = zeros(devnumber,1);
    AO.(ifam).Status      = zeros(devnumber,1);
    AO.(ifam).DeviceName  = cell(devnumber,1);
    AO.(ifam).CommonNames = cell(devnumber,1);
    AO.(ifam).Monitor.TangoNames  = cell(devnumber,1);
    % make Setpoint structure than specific data
    AO.(ifam).Setpoint = AO.(ifam).Monitor;
    AO.(ifam).Setpoint.Range = zeros(devnumber,2);
    for ik = 1: devnumber,
      AO.(ifam).ElementList(ik) = varlist.(ifam){ik,1}; 
      AO.(ifam).DeviceList(ik,:) = varlist.(ifam){ik,2};
      AO.(ifam).DeviceName(ik) = deblank(varlist.(ifam)(ik,3)); 
      AO.(ifam).Status(ik) = varlist.(ifam){ik,4};          
      AO.(ifam).CommonNames(ik) = deblank(varlist.(ifam)(ik,5));      
  
      AO.(ifam).Monitor.TangoNames(ik) = strcat(AO.(ifam).DeviceName{ik}, {'/currentPM'});

      AO.(ifam).Setpoint.TangoNames(ik) = strcat(AO.(ifam).DeviceName{ik},{'/currentPM'});
      AO.(ifam).Setpoint.Range(ik,:) = varlist.(ifam){ik,6};
    end

    
    %AO.(ifam).Setpoint = AO.(ifam).Monitor;
    
    AO.(ifam).Monitor.Handles(:,1) = NaN*ones(devnumber,1);
    
    % Conversion between HW and Physics
    HW2PhysicsParams = magnetcoefficients(AO.(ifam).FamilyName);  
     Physics2HWParams = magnetcoefficients(AO.(ifam).FamilyName); 
    val = 1.0; % scaling factor
    for ii=1:devnumber
      AO.(ifam).Monitor.HW2PhysicsParams{1}(ii,:) = HW2PhysicsParams';
      AO.(ifam).Monitor.HW2PhysicsParams{2}(ii,:) = val;
      AO.(ifam).Monitor.Physics2HWParams{1}(ii,:) = Physics2HWParams';
      AO.(ifam).Monitor.Physics2HWParams{2}(ii,:) = val;
    end 
    
    % same configuration for Monitor and Setpoint value concerning hardware to physics units
    AO.(ifam).Setpoint.HW2PhysicsParams = AO.(ifam).Monitor.HW2PhysicsParams;
    AO.(ifam).Setpoint.Physics2HWParams = AO.(ifam).Monitor.Physics2HWParams;
        
    AO.(ifam).Setpoint.MemberOf  = {'PlotFamily'};
    AO.(ifam).Setpoint.Tolerance     = 0.05;
        % information for getrunflag
    AO.(ifam).Setpoint.RunFlagFcn = @tangogetrunflag;
    AO.(ifam).Setpoint.RampRate = 1;
    
    AO.(ifam).Desired  = AO.(ifam).Monitor;
    
    
    %step size of k2 to correct the chromaticity
    %  used by stepchrom.m 
    deltak2 =1e-5;
    AO.(ifam).Setpoint.DeltaRespMat  = deltak2; % Physics units for a thick sextupole
    %need to reuse when the real machine is ready
    %AO.(ifam).Setpoint.DeltaRespMat=physics2hw(AO.(ifam).FamilyName,'Setpoint',AO.(ifam).Setpoint.DeltaRespMat,AO.(ifam).DeviceList);    
    
    
    %convert response matrix kicks to HWUnits (after AO is loaded to AppData)
    setao(AO);   %required to make physics2hw function
  

end
    
    % ifam = sprintf('SX%s', num2str(k));
    %devnumber = size(varlist.(ifam),1);
    %% preallocation
    % AO.(ifam).ElementList = zeros(devnumber,1);
    % AO.(ifam).Status      = zeros(devnumber,1);
    % AO.(ifam).DeviceName  = cell(devnumber,1);
    % AO.(ifam).DeviceName  = cell(devnumber,1);
    % AO.(ifam).CommonNames = cell(devnumber,1);
    % AO.(ifam).Monitor.TangoNames  = cell(devnumber,1);
    % % make Setpoint structure than specific data
    % AO.(ifam).Setpoint = AO.(ifam).Monitor;
    % AO.(ifam).Setpoint.Range = zeros(devnumber,2);
    
    % for ik = 1: devnumber,
    %     AO.(ifam).ElementList(ik)  = varlist.(ifam){ik,1};
    %     AO.(ifam).DeviceList(ik,:) = varlist.(ifam){ik,2};
    %     AO.(ifam).DeviceName(ik)   = deblank(varlist.(ifam)(ik,3));
    %     AO.(ifam).Status(ik)       = varlist.(ifam){ik,4};
    %     AO.(ifam).CommonNames(ik)  = deblank(varlist.(ifam)(ik,5));
    %     AO.(ifam).Monitor.TangoNames(ik)  = strcat(AO.(ifam).DeviceName{ik}, {'/currentPM'});
    %     AO.(ifam).Setpoint.TangoNames(ik) = strcat(AO.(ifam).DeviceName{ik}, {'/currentPM'});
    %     AO.(ifam).Setpoint.Range(ik,:)      = varlist.(ifam){ik,6};
    % end
    

%% All Sextupoles, configurations for the "plotfamily" feature
ifam = 'Sall';
AO.(ifam).FamilyName             = ifam;
AO.(ifam).MemberOf               = {ifam;'PlotFamily'};
AO.(ifam).FamilyType             = 'Sall';
AO.(ifam).Mode                   = Mode;

AO.(ifam).DeviceName= {};

for k = 1:3, 
    isext = sprintf('SX%d',k);
    AO.(ifam).DeviceName = {AO.(ifam).DeviceName{:} AO.(isext).DeviceName{:}};
end

AO.(ifam).DeviceName = AO.(ifam).DeviceName';

AO.(ifam).DeviceList = [];

devnumber = length(AO.(ifam).DeviceName);

% build fake device list in order to use plotfamily
for k =1:devnumber,
    AO.(ifam).DeviceList = [AO.(ifam).DeviceList; 1 k];
end

AO.(ifam).ElementList            = (1:devnumber)';
AO.(ifam).Status = ones(devnumber,1);
AO.(ifam).Monitor.Mode           = Mode;
AO.(ifam).Monitor.DataType       = 'Scalar';
AO.(ifam).Monitor.Units          = 'Hardware';
AO.(ifam).Monitor.HWUnits        = 'ampere';
AO.(ifam).Monitor.PhysicsUnits   = 'meter^-3';
AO.(ifam).Monitor.HW2PhysicsFcn  = @amp2k;
AO.(ifam).Monitor.Physics2HWFcn  = @k2amp;
AO.(ifam).Monitor.TangoNames     = strcat(AO.(ifam).DeviceName, '/current') ;    
AO.(ifam).Setpoint               = AO.(ifam).Monitor;
AO.(ifam).Setpoint.MemberOf      = {'PlotFamily'};
AO.(ifam).Setpoint.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/currentPM');

AO.(ifam).Voltage = AO.(ifam).Monitor;
AO.(ifam).Voltage.MemberOf  = {'PlotFamily'};
AO.(ifam).Voltage.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/voltage');

AO.(ifam).DCCT = AO.(ifam).Monitor;
AO.(ifam).DCCT.MemberOf  = {'PlotFamily'};
AO.(ifam).DCCT.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/current');

AO.(ifam).CurrentTotal = AO.(ifam).Monitor;
AO.(ifam).CurrentTotal.MemberOf  = {'PlotFamily'};
AO.(ifam).CurrentTotal.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/currentTotal');

AO.(ifam).CurrentSetpoint = AO.(ifam).Monitor;
AO.(ifam).CurrentSetpoint.MemberOf  = {'PlotFamily'};
AO.(ifam).CurrentSetpoint.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/currentSetpoint');

AO.(ifam).SumOffset = AO.(ifam).Monitor;
AO.(ifam).SumOffset.MemberOf  = {'PlotFamily'};
AO.(ifam).SumOffset.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/ecart3');

AO.(ifam).Offset1 = AO.(ifam).Monitor;
AO.(ifam).Offset1.MemberOf  = {'PlotFamily'};
AO.(ifam).Offset1.TangoNames(:,:)  = strcat(AO.(ifam).DeviceName,'/currentOffset1');

% Build up fake poistion for plotfamily
AO.(ifam).Position = (1:length(AO.(ifam).DeviceName))'*16.8/length(AO.(ifam).DeviceName);

%%%%%%%%%%%%%%%%%%
%% Pulsed Magnet
%%%%%%%%%%%%%%%%%%

%% All Injection kicker
ifam = 'K_INJ';
AO.(ifam).FamilyName           = ifam;
AO.(ifam).FamilyType           = ifam;
AO.(ifam).MemberOf             = {'Injection';'Archivable';'EP'};
AO.(ifam).Monitor.Mode         = Mode;
AO.(ifam).Monitor.DataType     = 'Scalar';

AO.(ifam).Status                   = ones(2,1);
AO.(ifam).DeviceName               = ['RI-C1/PE/KIC.010';'RI-C1/PE/KIC.020'];
AO.(ifam).CommonNames              = ifam;
AO.(ifam).ElementList              = [1 2]';
AO.(ifam).DeviceList(:,:)          = [1 1; 1 2];
AO.(ifam).Monitor.Handles(:,1)     = NaN*ones(2,1);
AO.(ifam).Monitor.TangoNames       = strcat(AO.(ifam).DeviceName, '/voltage');
AO.(ifam).Monitor.Units            = 'Hardware';
AO.(ifam).Monitor.HWUnits          = 'V';
AO.(ifam).Monitor.PhysicsUnits     = 'mrad';
AO.(ifam).Monitor.HW2PhysicsFcn     = @amp2k;
AO.(ifam).Monitor.Physics2HWFcn     = @k2amp;
HW2PhysicsParams                    = magnetcoefficients(ifam);
Physics2HWParams                    = HW2PhysicsParams;
AO.(ifam).Monitor.HW2PhysicsParams{1}(:,:) = HW2PhysicsParams;
AO.(ifam).Monitor.HW2PhysicsParams{2}(:,:) = 1;
AO.(ifam).Monitor.Physics2HWParams{1}(:,:) = Physics2HWParams;
AO.(ifam).Monitor.Physics2HWParams{2}(:,:) = 1;
AO.(ifam).Setpoint = AO.(ifam).Monitor;
AO.(ifam).Desired = AO.(ifam).Monitor;



%% Septum 
ifam = 'SEP';
AO.(ifam).FamilyName           = ifam;
AO.(ifam).FamilyType           = ifam;
AO.(ifam).MemberOf             = {'Injection';'Archivable';'EP'};
AO.(ifam).Monitor.Mode         = Mode;
AO.(ifam).Monitor.DataType     = 'Scalar';

AO.(ifam).Status                   = 1;
AO.(ifam).DeviceName               = cellstr('RI-C1/ME/SEP.010');
AO.(ifam).CommonNames              = ifam;
AO.(ifam).ElementList              = 1;
AO.(ifam).DeviceList(:,:)           = [1 1];
AO.(ifam).Monitor.Handles(:,1)     = NaN;
AO.(ifam).Monitor.TangoNames       = strcat(AO.(ifam).DeviceName, '/voltage');
AO.(ifam).Monitor.HW2PhysicsFcn     = @amp2k;
AO.(ifam).Monitor.Physics2HWFcn     = @k2amp;
HW2PhysicsParams                    = magnetcoefficients(ifam);
Physics2HWParams                    = HW2PhysicsParams;
AO.(ifam).Monitor.HW2PhysicsParams{1}(:,:) = HW2PhysicsParams;
AO.(ifam).Monitor.HW2PhysicsParams{2}(:,:) = 1;
AO.(ifam).Monitor.Physics2HWParams{1}(:,:) = Physics2HWParams;
AO.(ifam).Monitor.Physics2HWParams{2}(:,:) = 1;
AO.(ifam).Monitor.Units            = 'Hardware';
AO.(ifam).Monitor.HWUnits          = 'V';
AO.(ifam).Monitor.PhysicsUnits     = 'mrad';
AO.(ifam).Setpoint = AO.(ifam).Monitor;
AO.(ifam).Desired = AO.(ifam).Monitor;


%%%%%%%%%%%%%%%%%%%%
%%tune correctors
%%%%%%%%%%%%%%%%%%%%
AO.QP31.MemberOf = {AO.QP31.MemberOf{:}  'Tune Corrector'}';
AO.QP4.MemberOf = {AO.QP4.MemberOf{:}  'Tune Corrector'}';

%% chromaticity correctors
AO.SX1.MemberOf  = {AO.SX1.MemberOf{:} 'Chromaticity Corrector'}';
AO.SX2.MemberOf = {AO.SX2.MemberOf{:} 'Chromaticity Corrector'}';
%%%%%%%%%%%%%%%%%%%%

%% CYCLAGE
%%%%%%%%%%%%%%%%%
disp('cycling configuration ...');

%% cycleramp For dipole magnet
ifam = 'CycleBEND';

AO.(ifam).FamilyName             = 'CycleBEND';
AO.(ifam).MemberOf               = {'CycleBEND'; 'Cyclage'};
AO.(ifam).Mode                   = Mode;


AO.(ifam).DeviceName             = 'ANS/AE/cycleDipole';
AO.(ifam).DeviceList             = [1 1];
AO.(ifam).ElementList            = 1;
AO.(ifam).Inom = 541.789;
AO.(ifam).Imax = 579.9;
AO.(ifam).Status = 1;

if ControlRoomFlag
        AO.(ifam).GroupId                = tango_group_create(ifam);
        tango_group_add(AO.(ifam).GroupId, AO.(ifam).DeviceName);    
else
    AO.(ifam).GroupId = nan;
end

%% cycleramp For H-corrector magnets
ifam  = 'CycleHCOR';
ifamQ = 'HCOR';

AO.(ifam).FamilyName             = ifam;
AO.(ifam).MemberOf               = {'CycleCOR'; ifam; 'Cyclage'};
AO.(ifam).Mode                   = Mode;
dev = getfamilydata(ifamQ,'DeviceName');
AO.(ifam).DeviceName             = regexprep(dev,'/S','/cycleS');
AO.(ifam).DeviceList             = AO.(ifamQ).DeviceList;
AO.(ifam).ElementList            = AO.(ifamQ).ElementList;

if ControlRoomFlag
        AO.(ifam).GroupId                = tango_group_create(ifam);
        tango_group_add(AO.(ifam).GroupId, AO.(ifam).DeviceName(find(AO.(ifamQ).Status),:)')
else
    AO.(ifam).GroupId = nan;
end

devnumber = length(AO.(ifam).DeviceName);
AO.(ifam).Inom = 1.*ones(1,devnumber);
AO.(ifam).Imax = 10.99*ones(1,devnumber);
AO.(ifam).Status = AO.(ifamQ).Status;
AO.(ifam).Monitor.Mode           = Mode;
AO.(ifam).Monitor.Handles(:,1)   = NaN*ones(devnumber,1);
AO.(ifam).Monitor.DataType       = 'Scalar';
AO.(ifam).Monitor.Units          = 'Hardware';
AO.(ifam).Monitor.HWUnits        = 'ampere';
AO.(ifam).Monitor.PhysicsUnits   = 'rad';
AO.(ifam).Monitor.TangoNames   = strcat(AO.(ifam).DeviceName, '/totalProgression');

%% cycleramp For V-corrector magnets
ifam = 'CycleVCOR';
ifamQ = 'VCOR';

AO.(ifam).FamilyName             = ifam;
AO.(ifam).MemberOf               = {'CycleCOR'; ifam; 'Cyclage'};
AO.(ifam).Mode                   = Mode;
dev = getfamilydata(ifamQ,'DeviceName');
AO.(ifam).DeviceName             = regexprep(dev,'/S','/cycleS');
AO.(ifam).DeviceList             = AO.(ifamQ).DeviceList;
AO.(ifam).ElementList            = AO.(ifamQ).ElementList;


if ControlRoomFlag
    AO.(ifam).GroupId                = tango_group_create(ifam);
    tango_group_add(AO.(ifam).GroupId, AO.(ifam).DeviceName(find(AO.(ifamQ).Status),:)')
else
    AO.(ifam).GroupId = nan;
end

devnumber = length(AO.(ifam).DeviceName);
AO.(ifam).Inom = 1.*ones(1,devnumber);
AO.(ifam).Imax = 13.99*ones(1,devnumber);
AO.(ifam).Status = AO.(ifamQ).Status;
AO.(ifam).Monitor.Mode           = Mode;
AO.(ifam).Monitor.Handles(:,1)   = NaN*ones(devnumber,1);
AO.(ifam).Monitor.DataType       = 'Scalar';
AO.(ifam).Monitor.Units          = 'Hardware';
AO.(ifam).Monitor.HWUnits        = 'ampere';
AO.(ifam).Monitor.PhysicsUnits   = 'rad';
AO.(ifam).Monitor.TangoNames   = strcat(AO.(ifam).DeviceName, '/totalProgression');

% cycleramp For quadrupoles magnets
%% CYCLEQP1
ifam = 'CycleQP1';
ifamQ = 'QP1';

AO.(ifam).FamilyName             = ifam;
AO.(ifam).MemberOf               = {ifam; 'CycleQP'; 'Cyclage'};
AO.(ifam).Mode                   = Mode;
dev = getfamilydata(ifamQ,'DeviceName');
AO.(ifam).DeviceName             = regexprep(dev,'/QP1','/cycleQP');
AO.(ifam).DeviceList             = family2dev(ifamQ);
AO.(ifam).ElementList            = family2elem(ifamQ);

if ControlRoomFlag
        %add devices to group
        AO.(ifam).GroupId                = tango_group_create(ifam);
        tango_group_add(AO.(ifam).GroupId, AO.(ifam).DeviceName');
else
    AO.(ifam).GroupId = nan;
end

devnumber = length(AO.(ifam).ElementList);
%AO.(ifam).Inom = 200.*ones(1,devnumber);
AO.(ifam).Imax = -249*ones(1,devnumber);
AO.(ifam).Status = ones(devnumber,1);
AO.(ifam).Monitor.Mode           = Mode;
AO.(ifam).Monitor.Handles(:,1)   = NaN*ones(devnumber,1);
AO.(ifam).Monitor.DataType       = 'Scalar';
AO.(ifam).Monitor.Units          = 'Hardware';
AO.(ifam).Monitor.HWUnits        = 'ampere';
AO.(ifam).Monitor.PhysicsUnits   = 'rad';
AO.(ifam).Monitor.TangoNames   = strcat(AO.(ifam).DeviceName, '/totalProgression');

%% CYCLEQ2
ifam = 'CycleQP2';
ifamQ = 'QP2';

AO.(ifam).FamilyName             = ifam;
AO.(ifam).MemberOf               = {ifam;'CycleQP2'; 'Cyclage'};
AO.(ifam).Mode                   = Mode;
dev = getfamilydata(ifamQ,'DeviceName');
AO.(ifam).DeviceName             = regexprep(dev,'/QP','/cycleQP');
AO.(ifam).DeviceList             = family2dev(ifamQ);
AO.(ifam).ElementList            = family2elem(ifamQ);

if ControlRoomFlag
        %add devices to group
        AO.(ifam).GroupId                = tango_group_create(ifam);
        tango_group_add(AO.(ifam).GroupId, AO.(ifam).DeviceName');
else
    AO.(ifam).GroupId = nan;
end

devnumber = length(AO.(ifam).ElementList);
%AO.(ifam).Inom = 150*ones(1,devnumber);
AO.(ifam).Imax = 249*ones(1,devnumber);
AO.(ifam).Status = ones(devnumber,1);
AO.(ifam).Monitor.Mode           = Mode;
AO.(ifam).Monitor.Handles(:,1)   = NaN*ones(devnumber,1);
AO.(ifam).Monitor.DataType       = 'Scalar';
AO.(ifam).Monitor.Units          = 'Hardware';
AO.(ifam).Monitor.HWUnits        = 'ampere';
AO.(ifam).Monitor.PhysicsUnits   = 'rad';
AO.(ifam).Monitor.TangoNames   = strcat(AO.(ifam).DeviceName, '/totalProgression');

%% CYCLEQP3
ifam = 'CycleQP3';
ifamQ = 'QP3';

AO.(ifam).FamilyName             = ifam;
AO.(ifam).MemberOf               = {ifam;'CycleQP3'; 'Cyclage'};
AO.(ifam).Mode                   = Mode;
dev = getfamilydata(ifamQ,'DeviceName');
AO.(ifam).DeviceName             = regexprep(dev,'/QP','/cycleQP');
AO.(ifam).DeviceList             = family2dev(ifamQ);
AO.(ifam).ElementList            = family2elem(ifamQ);

if ControlRoomFlag
        %add devices to group
        AO.(ifam).GroupId                = tango_group_create(ifam);
        tango_group_add(AO.(ifam).GroupId, AO.(ifam).DeviceName');
else
    AO.(ifam).GroupId = nan;
end

devnumber = length(AO.(ifam).ElementList);
%AO.(ifam).Inom = 150*ones(1,devnumber);
AO.(ifam).Imax = -249*ones(1,devnumber);
AO.(ifam).Status = ones(devnumber,1);
AO.(ifam).Monitor.Mode           = Mode;
AO.(ifam).Monitor.Handles(:,1)   = NaN*ones(devnumber,1);
AO.(ifam).Monitor.DataType       = 'Scalar';
AO.(ifam).Monitor.Units          = 'Hardware';
AO.(ifam).Monitor.HWUnits        = 'ampere';
AO.(ifam).Monitor.PhysicsUnits   = 'rad';
AO.(ifam).Monitor.TangoNames   = strcat(AO.(ifam).DeviceName, '/totalProgression');

%% CYCLEQP4
ifam = 'CycleQP4';
ifamQ = 'QP4';

AO.(ifam).FamilyName             = ifam;
AO.(ifam).MemberOf               = {ifam;'CycleQP4'; 'Cyclage'};
AO.(ifam).Mode                   = Mode;
dev = getfamilydata(ifamQ,'DeviceName');
AO.(ifam).DeviceName             = regexprep(dev,'/QP','/cycleQP');
AO.(ifam).DeviceList             = family2dev(ifamQ);
AO.(ifam).ElementList            = family2elem(ifamQ);
if ControlRoomFlag
        %add devices to group
        AO.(ifam).GroupId                = tango_group_create(ifam);
        tango_group_add(AO.(ifam).GroupId, AO.(ifam).DeviceName');
else
    AO.(ifam).GroupId = nan;
end
devnumber = length(AO.(ifam).ElementList);
%AO.(ifam).Inom = 150*ones(1,devnumber);
AO.(ifam).Imax = -249*ones(1,devnumber);
AO.(ifam).Status = ones(devnumber,1);
AO.(ifam).Monitor.Mode           = Mode;
AO.(ifam).Monitor.Handles(:,1)   = NaN*ones(devnumber,1);
AO.(ifam).Monitor.DataType       = 'Scalar';
AO.(ifam).Monitor.Units          = 'Hardware';
AO.(ifam).Monitor.HWUnits        = 'ampere';
AO.(ifam).Monitor.PhysicsUnits   = 'rad';
AO.(ifam).Monitor.TangoNames   = strcat(AO.(ifam).DeviceName, '/totalProgression');

%% CYCLEQP31
ifam = 'CycleQP31';
ifamQ = 'QP31';

AO.(ifam).FamilyName             = ifam;
AO.(ifam).MemberOf               = {ifam;'CycleQP31'; 'Cyclage'};
AO.(ifam).Mode                   = Mode;
dev = getfamilydata(ifamQ,'DeviceName');
AO.(ifam).DeviceName             = regexprep(dev,'/QP','/cycleQP');
AO.(ifam).DeviceList             = family2dev(ifamQ);
AO.(ifam).ElementList            = family2elem(ifamQ);

if ControlRoomFlag
        %add devices to group
        AO.(ifam).GroupId                = tango_group_create(ifam);
        tango_group_add(AO.(ifam).GroupId, AO.(ifam).DeviceName');
else
    AO.(ifam).GroupId = nan;
end
devnumber = length(AO.(ifam).ElementList);
%AO.(ifam).Inom = 150*ones(1,devnumber);
AO.(ifam).Imax = 249*ones(1,devnumber);
AO.(ifam).Status = ones(devnumber,1);
AO.(ifam).Monitor.Mode           = Mode;
AO.(ifam).Monitor.Handles(:,1)   = NaN*ones(devnumber,1);
AO.(ifam).Monitor.DataType       = 'Scalar';
AO.(ifam).Monitor.Units          = 'Hardware';
AO.(ifam).Monitor.HWUnits        = 'ampere';
AO.(ifam).Monitor.PhysicsUnits   = 'rad';
AO.(ifam).Monitor.TangoNames   = strcat(AO.(ifam).DeviceName, '/totalProgression');

%% CYCLEQP41
ifam = 'CycleQP41';
ifamQ = 'QP41';

AO.(ifam).FamilyName             = ifam;
AO.(ifam).MemberOf               = {ifam;'CycleQP41'; 'Cyclage'};
AO.(ifam).Mode                   = Mode;
dev = getfamilydata(ifamQ,'DeviceName');
AO.(ifam).DeviceName             = regexprep(dev,'/QP','/cycleQP');
AO.(ifam).DeviceList             = family2dev(ifamQ);
AO.(ifam).ElementList            = family2elem(ifamQ);
if ControlRoomFlag
        %add devices to group
        AO.(ifam).GroupId                = tango_group_create(ifam);
        tango_group_add(AO.(ifam).GroupId, AO.(ifam).DeviceName');
else
    AO.(ifam).GroupId = nan;
end
devnumber = length(AO.(ifam).ElementList);
%AO.(ifam).Inom = 150*ones(1,devnumber);
AO.(ifam).Imax = -249*ones(1,devnumber);
AO.(ifam).Status = ones(devnumber,1);
AO.(ifam).Monitor.Mode           = Mode;
AO.(ifam).Monitor.Handles(:,1)   = NaN*ones(devnumber,1);
AO.(ifam).Monitor.DataType       = 'Scalar';
AO.(ifam).Monitor.Units          = 'Hardware';
AO.(ifam).Monitor.HWUnits        = 'ampere';
AO.(ifam).Monitor.PhysicsUnits   = 'rad';
AO.(ifam).Monitor.TangoNames   = strcat(AO.(ifam).DeviceName, '/totalProgression');

%% cycleramp for sextupole magnets
%% CYCLES1
ifam = 'CycleSX1';
ifamQ = 'SX1';

AO.(ifam).FamilyName             = ifam;
AO.(ifam).MemberOf               = {ifam; 'Cyclage'};
AO.(ifam).Mode                   = Mode;
dev = getfamilydata(ifamQ,'DeviceName');
AO.(ifam).DeviceName             = regexprep(dev,'/SX','/cycleSX');
AO.(ifam).DeviceList             = family2dev(ifamQ);
AO.(ifam).ElementList            = family2elem(ifamQ);
if ControlRoomFlag
        %add devices to group
        AO.(ifam).GroupId        = tango_group_create(ifam);
        tango_group_add(AO.(ifam).GroupId, AO.(ifam).DeviceName');
else
    AO.(ifam).GroupId = nan;
end
devnumber = length(AO.(ifam).ElementList);
AO.(ifam).Imax = 349*ones(1,devnumber);
AO.(ifam).Status = ones(devnumber,1);
AO.(ifam).Monitor.Mode           = Mode;
AO.(ifam).Monitor.Handles(:,1)   = NaN*ones(devnumber,1);
AO.(ifam).Monitor.DataType       = 'Scalar';
AO.(ifam).Monitor.Units          = 'Hardware';
AO.(ifam).Monitor.HWUnits        = 'ampere';
AO.(ifam).Monitor.PhysicsUnits   = 'rad';
AO.(ifam).Monitor.TangoNames   = strcat(AO.(ifam).DeviceName, '/totalProgression');

%% CYCLES2
ifam = 'CycleSX2';
ifamQ = 'SX2';

AO.(ifam).FamilyName             = ifam;
AO.(ifam).MemberOf               = {ifam; 'Cyclage'};
AO.(ifam).Mode                   = Mode;
dev = getfamilydata(ifamQ,'DeviceName');
AO.(ifam).DeviceName             = regexprep(dev,'/SX','/cycleSX');
AO.(ifam).DeviceList             = family2dev(ifamQ);
AO.(ifam).ElementList            = family2elem(ifamQ);
if ControlRoomFlag
        %add devices to group
        AO.(ifam).GroupId        = tango_group_create(ifam);
        tango_group_add(AO.(ifam).GroupId, AO.(ifam).DeviceName');
else
    AO.(ifam).GroupId = nan;
end
devnumber = length(AO.(ifam).ElementList);
AO.(ifam).Imax = -349*ones(1,devnumber);
AO.(ifam).Status = ones(devnumber,1);
AO.(ifam).Monitor.Mode           = Mode;
AO.(ifam).Monitor.Handles(:,1)   = NaN*ones(devnumber,1);
AO.(ifam).Monitor.DataType       = 'Scalar';
AO.(ifam).Monitor.Units          = 'Hardware';
AO.(ifam).Monitor.HWUnits        = 'ampere';
AO.(ifam).Monitor.PhysicsUnits   = 'rad';
AO.(ifam).Monitor.TangoNames   = strcat(AO.(ifam).DeviceName, '/totalProgression');

%% CYCLES3
ifam = 'CycleSX3';
ifamQ = 'SX3';

AO.(ifam).FamilyName             = ifam;
AO.(ifam).MemberOf               = {ifam; 'Cyclage'};
AO.(ifam).Mode                   = Mode;
dev = getfamilydata(ifamQ,'DeviceName');
AO.(ifam).DeviceName             = regexprep(dev,'/SX','/cycleSX');
AO.(ifam).DeviceList             = family2dev(ifamQ);
AO.(ifam).ElementList            = family2elem(ifamQ);
if ControlRoomFlag
        %add devices to group
        AO.(ifam).GroupId        = tango_group_create(ifam);
        tango_group_add(AO.(ifam).GroupId, AO.(ifam).DeviceName');
else
    AO.(ifam).GroupId = nan;
end
devnumber = length(AO.(ifam).ElementList);
AO.(ifam).Imax = -349*ones(1,devnumber);
AO.(ifam).Status = ones(devnumber,1);
AO.(ifam).Monitor.Mode           = Mode;
AO.(ifam).Monitor.Handles(:,1)   = NaN*ones(devnumber,1);
AO.(ifam).Monitor.DataType       = 'Scalar';
AO.(ifam).Monitor.Units          = 'Hardware';
AO.(ifam).Monitor.HWUnits        = 'ampere';
AO.(ifam).Monitor.PhysicsUnits   = 'rad';
AO.(ifam).Monitor.TangoNames   = strcat(AO.(ifam).DeviceName, '/totalProgression');

%============
%% RF System
%============
ifam = 'RF';
AO.(ifam).FamilyName                = ifam;
AO.(ifam).FamilyType                = 'RF';
AO.(ifam).MemberOf                  = {'RF','RFSystem'};
AO.(ifam).Status                    = 1;
AO.(ifam).CommonNames               = 'RF';
AO.(ifam).DeviceList                = [1 1];
AO.(ifam).ElementList               = 1;
AO.(ifam).DeviceName(:,:)           = {'RI-C2/RF/RFC.010/'};

%Frequency Readback
AO.(ifam).Monitor.Mode                = Mode;
% AO.(ifam).Monitor.Mode                = 'Special';
% AO.(ifam).Monitor.SpecialFunctionGet = 'getrf2';
AO.(ifam).Monitor.DataType            = 'Scalar';
AO.(ifam).Monitor.Units               = 'Hardware';
AO.(ifam).Monitor.HW2PhysicsParams    = 1e+6;       %no hw2physics function necessary
AO.(ifam).Monitor.Physics2HWParams    = 1e-6;
AO.(ifam).Monitor.HWUnits             = 'MHz';
AO.(ifam).Monitor.PhysicsUnits        = 'Hz';
AO.(ifam).Monitor.TangoNames          = {'ANS/RF/MasterClock/frequency'};
AO.(ifam).Monitor.Handles             = NaN;
AO.(ifam).Monitor.Range               = [400 600];

AO.(ifam).Setpoint = AO.(ifam).Monitor;
AO.(ifam).Desired  = AO.(ifam).Monitor;

%============
%% Cryomodules System
%=========

%=======
%% CHROMATICITIES: Soft Family
%=======
ifam = 'CHRO';
AO.(ifam).FamilyName  = ifam;
AO.(ifam).FamilyType  = 'Diagnostic';
AO.(ifam).MemberOf    = {'Diagnostics', 'CHRO'};
AO.(ifam).CommonNames = ['xix';'xiz'];
AO.(ifam).DeviceList  = [ 1 1; 1 2;];
AO.(ifam).ElementList = [1; 2];
AO.(ifam).Status      = [1; 1];
%======
%% Dipole B-FIELD from RMN probe
%======

%=======
%% TUNE
%=======
ifam = 'TUNE';
AO.(ifam).FamilyName  = ifam;
AO.(ifam).FamilyType  = 'Diagnostic';
AO.(ifam).MemberOf    = {'Diagnostics', 'TUNE'};
AO.(ifam).CommonNames = ['nux';'nuz';'nus'];
AO.(ifam).DeviceList  = [1 1; 1 2; 1 3];
AO.(ifam).ElementList = [1 2 3]';
AO.(ifam).Status      = [1 1 1]';

AO.(ifam).Monitor.Mode                   = Mode; %'Simulator';  % Mode;
AO.(ifam).Monitor.DataType               = 'Scalar';
AO.(ifam).Monitor.DataTypeIndex          = [1 2];

% need to customized for THOMX
AO.(ifam).Monitor.TangoNames             = ['ANS/DG/BPM-TUNEX/Nu';'ANS/DG/BPM-TUNEZ/Nu';'ANS/DG/BPM-TUNEZ/Nu'];

AO.(ifam).Monitor.Units                  = 'Hardware';
AO.(ifam).Monitor.Handles                = NaN;
AO.(ifam).Monitor.HW2PhysicsParams       = 1;
AO.(ifam).Monitor.Physics2HWParams       = 1;
AO.(ifam).Monitor.HWUnits                = 'fractional tune';
AO.(ifam).Monitor.PhysicsUnits           = 'fractional tune';
%=======
%% TUNEFBT
% control
%=======

%=======
%% Coupling
%=======

%======
%% DCCT
% measurement (Tango)
%======

% Need to check in the future...
ifam = 'DCCT';
AO.(ifam).FamilyName                     = ifam;
AO.(ifam).FamilyType                     = 'Diagnostic';
AO.(ifam).MemberOf                       = {'Diagnostics','DCCT'};
AO.(ifam).CommonNames                    = 'DCCT';
AO.(ifam).DeviceList                     = [1 1];
AO.(ifam).ElementList                    = 1;
AO.(ifam).Status                         = AO.(ifam).ElementList;

AO.(ifam).Monitor.Mode                   = Mode;
AO.(ifam).FamilyName                     = 'DCCT';
AO.(ifam).Monitor.DataType               = 'Scalar';
AO.(ifam).Monitor.TangoNames             = 'ANS/DG/DCCT-CTRL/current'; %afin de ne pas avoir de bug
AO.(ifam).Monitor.Units                  = 'Hardware';
AO.(ifam).Monitor.Handles                = NaN;
AO.(ifam).Monitor.HWUnits                = 'mA';
AO.(ifam).Monitor.PhysicsUnits           = 'A';
AO.(ifam).Monitor.HW2PhysicsParams       = 1;
AO.(ifam).Monitor.Physics2HWParams       = 1;


%======
%% Lifetime
% measurement
%======
    ifam = 'LifeTime';
AO.(ifam).FamilyName                     = ifam;
AO.(ifam).FamilyType                     = 'Diagnostic';
AO.(ifam).MemberOf                       = {'Diagnostics','DCCT'};
AO.(ifam).CommonNames                    = 'LifeTime';
AO.(ifam).DeviceList                     = [1 1];
AO.(ifam).ElementList                    = 1;
AO.(ifam).Status                         = AO.(ifam).ElementList;

AO.(ifam).Monitor.Mode                   = Mode;
AO.(ifam).FamilyName                     = 'LifeTime';
AO.(ifam).Monitor.DataType               = 'Scalar';
AO.(ifam).Monitor.TangoNames             = 'ANS/DG/DCCT-CTRL/lifeTime'; %afin de ne pas avoir de bug
AO.(ifam).Monitor.Units                  = 'Hardware';
AO.(ifam).Monitor.Handles                = NaN;
AO.(ifam).Monitor.HWUnits                = 'Seconds';
AO.(ifam).Monitor.PhysicsUnits           = 'Seconds';
AO.(ifam).Monitor.HW2PhysicsParams       = 1;
AO.(ifam).Monitor.Physics2HWParams       = 1;

%%%%%%%%%%%%%%%%%%
%% Alignment error 
    % setting
%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%
%% AVLS GEOPHONES
%  ground vibration setting
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AVLS GEOPHONES PEAK VALUE
%  ground vibration setting
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
%% VACUUM SYSTEM
%%%%%%%%%%%%%%%%%%

%=====================
%% Gamma Monitors DOSE
% BLM set & read
%====================

%=========================
%% Gamma Monitors DOSErate
% BLM set & read
%=========================
    
%=====================
%% Neutron Monitors DOSE
%====================

%============================
%% Neutron Monitors DOSE RATE
%============================

% %====================
% %% Machine Parameters
% %====================

%======
%% Septum
%======


% Save AO
setao(AO);

% The operational mode sets the path, filenames, and other important params
% Run setoperationalmode after most of the AO is built so that the Units and Mode fields
% can be set in setoperationalmode

waitbar(0.4,h);

setoperationalmode(OperationalMode);
%run('/home/operateur/GrpPhysiqueMachine/Laurent/matlab/nano/setdevelopmentmode.m')
%run([getdvptdirectory '/setdevelopmentmode.m']);
%run('/home/operateur/temp/testMML/setdevelopmentmode.m');

waitbar(0.5,h);

%======================================================================
%======================================================================
%% Append Accelerator Toolbox information
%======================================================================
%======================================================================
disp('** Initializing Accelerator Toolbox information');

AO = getao;

%% Machine Params
ifam = ('MachineParams');
AO.(ifam).AT.ATType       = 'MachineParams';
AO.(ifam).AT.ATName(1,:)  = 'Energy  ';
AO.(ifam).AT.ATName(2,:)  = 'current ';
AO.(ifam).AT.ATName(3,:)  = 'Lifetime';

% Save AO
setao(AO);

disp('Setting min max configuration from TANGO static database ...');

waitbar(0.80,h);

disp('Setting gain offset configuration  ...');

setfamilydata([0.175; 0.720; NaN],'TUNE','Golden');
setfamilydata([0.0; 0.0],'CHRO','Golden');
setfamilydata_local('BPMx');
setfamilydata_local('BPMz');
setfamilydata_local('HCOR');
setfamilydata_local('VCOR');
setfamilydata_local('QP1');
setfamilydata_local('QP2');
setfamilydata_local('QP31');
setfamilydata_local('QP41');
setfamilydata_local('QP4');
setfamilydata_local('QP3');

setfamilydata_local('SX1');
setfamilydata_local('SX2');
setfamilydata_local('SX3');

waitbar(0.95,h);

if iscontrolroom
    switch2online;
else
    switch2sim;
end

delete(h);
% set close orbit
% fprintf('%3d %3d  %10.6f  %10.6f\n', [family2dev('BPMx'), getgolden('BPMx'), getgolden('BPMz')]');
% fprintf('%3d %3d  %10.6f  %10.6f\n', [family2dev('BPMx'), getam('BPMx'), getam('BPMz')]');
% fprintf('%3d %3d  %10.6f  %10.6f\n', [family2dev('BPMx'), gethbpmaverage, getvbpmaverage]'); % With averaging




  %% 20 janvier 2012, 9mA, 122 BPM, 1/4 RUN5
% Golden = [
%setfamilydata(Golden(:,3),'BPMx','Golden',Golden(:,1:2));
%setfamilydata(Golden(:,4),'BPMz','Golden',Golden(:,1:2));

end
  
%% LOCAL FUNCTIONS
function local_tango_kill_allgroup(AO)
% kill all group if exist
FamilyList = getfamilylist('Cell');
for k=1:length(FamilyList)
    if isfield(AO.(FamilyList{k}), 'GroupId'),
        tango_group_kill(AO.(FamilyList{k}).GroupId);
    end
end
end

function setfamilydata_local(Family)
% set all data in one command

if ismemberof(Family,'QUAD') || ismemberof(Family,'SEXT') || ...
         ismemberof(Family,'BEND')
    setfamilydata(1,Family,'Gain');
    setfamilydata(0,Family,'Offset');
    setfamilydata(0,Family,'Coupling');
end

if ismemberof(Family,'BPM')
    setfamilydata(0.001,Family,'Sigma');
    setfamilydata(0.0,Family,'Golden');   
    %setfamilydata(measdisp(Family,'struct','model'),Family,'Dispersion'); % needed for orbit correction w/ RF
end

% Order fields for all families
AO = getao;
Familylist = getfamilylist;
for k = 1:length(Familylist),
    FamilyName = deblank(Familylist(k,:));
    AO.(FamilyName) = orderfields(AO.(FamilyName));
    FieldNamelist = fieldnames(AO.(FamilyName));
    for k1 = 1:length(FieldNamelist);
        FieldName = FieldNamelist{k1};
        if isstruct((AO.(FamilyName).(FieldName)))
            AO.(FamilyName).(FieldName) = orderfields(AO.(FamilyName).(FieldName));
        end
    end
end
setao(AO)
end