function boosterinit(OperationalMode)
%BOOSTERINIT - Initializes parameters for SOLEIL Booster control in MATLAB
%
% Written by Laurent S. Nadolski, Synchrotron SOLEIL
% 
%==========================
% Accelerator Family Fields
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
% HW2PhysicsParams      parameters used for conversion function
% Physics2HWFcn         function handle used to convert from physics to hardware units
% Physics2HWParams      parameters used for conversion function
% HWUnits               units for Hardware 'A';
% PhysicsUnits          units for physics 'Rad';
% Handles               monitor handle
%
% SETPOINT FIELDS
% Mode                  online/manual/special/simulator
% TangoNames            Devices tango names
% Units                 hardware or physics
% HW2PhysicsFcn         function handle used to convert from hardware to physics units
% HW2PhysicsParams      parameters used for conversion function
% Physics2HWFcn         function handle used to convert from physics to hardware units
% Physics2HWParams      parameters used for conversion function
% HWUnits               units for Hardware 'A';
% PhysicsUnits          units for physics 'Rad';
% Range                 minsetpoint, maxsetpoint;
% Tolerance             setpoint-monitor
% Handles               setpoint handle
%
%=============================================
% Accelerator Toolbox Simulation Fields
%=============================================
% ATType                Quad, Sext, etc
% ATIndex               index in THERING
% ATParameterGroup      parameter group
%
%============
% Family List
%============
%    BPMx
%    BPMz
%    HCOR
%    VCOR
%    BEND
%    QF and QD
%    SF and SD
%    RF
%    TUNE
%    DCCT
%    Machine Parameters
% If controlromm user is operator and online mode

[statuss WHO] = system('whoami');
% system gives back an visible character: carriage return!
% so comparison on the number of caracters
if strncmp(WHO, 'operateur',9),
    ControlRoomFlag = 1;
    Mode = 'Online';
else
    ControlRoomFlag = 0;
    Mode = 'Simulator';
end
 
%==============================
%load AcceleratorData structure
%==============================

if nargin < 1
    OperationalMode = 1;
end

global GLOBVAL THERING


Mode             = 'Online';
setad([]);       %clear AcceleratorData memory
AD.SubMachine       = 'Booster'; % Machine Name
AD.Energy        = 2.75; % Energy in GeV


setad(AD);   %load AcceleratorData

%%%%%%%%%%%%%%%%%%%%
% ACCELERATOR OBJECT
%%%%%%%%%%%%%%%%%%%%

setao([]);   %clear previous AcceleratorObjects

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BPMx data: status field designates if BPM in use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AO{1}.FamilyName               = 'BPMx';
AO{1}.FamilyType               = 'BPM';
AO{1}.MemberOf                 = {'BPM'; 'HBPM'; 'Diagnostics'};
AO{1}.Monitor.Mode             = Mode;
AO{1}.Monitor.Units            = 'Hardware';
AO{1}.Monitor.HWUnits          = 'mm';
AO{1}.Monitor.PhysicsUnits     = 'm';

% Get mapping from TANGO static database
if ControlRoomFlag
    map = tango_get_db_property('booster','tracy_bpm_mapping');
else
    map={'BPM01::BOO-C01/DG/BPM.01'
        'BPM02::BOO-C03/DG/BPM.02'
        'BPM03::BOO-C03/DG/BPM.03'
        'BPM04::BOO-C05/DG/BPM.04'
        'BPM05::BOO-C05/DG/BPM.05'
        'BPM06::BOO-C06/DG/BPM.06'
        'BPM07::BOO-C07/DG/BPM.07'
        'BPM08::BOO-C08/DG/BPM.08'
        'BPM09::BOO-C09/DG/BPM.09'
        'BPM10::BOO-C10/DG/BPM.10'
        'BPM11::BOO-C11/DG/BPM.11'
        'BPM12::BOO-C12/DG/BPM.12'
        'BPM13::BOO-C14/DG/BPM.13'
        'BPM14::BOO-C14/DG/BPM.14'
        'BPM15::BOO-C16/DG/BPM.15'
        'BPM16::BOO-C16/DG/BPM.16'
        'BPM17::BOO-C17/DG/BPM.17'
        'BPM18::BOO-C18/DG/BPM.18'
        'BPM19::BOO-C19/DG/BPM.19'
        'BPM20::BOO-C20/DG/BPM.20'
        'BPM21::BOO-C21/DG/BPM.21'
        'BPM22::BOO-C22/DG/BPM.22'};
    map=map';
end

cellindex = cell2mat(regexpi(map,'C[0-9]','once'))+1;
numindex  = cell2mat(regexpi(map,'\.[0-9]','once'))+1;
elemindex = cell2mat(regexpi(map,'BPM[0-9]','once'))+3;
sep       = cell2mat(regexpi(map,'::','once'))-1;
% dev       = deblank(regexprep(char(map),'^BPM\d*::',''));
dev       = regexprep(map,'^BPM\d*::','')';

nb = size(map,2);
ik = 1;
prev0 = 0;
for k = 1:nb
    prev = str2num(map{k}(cellindex(k):cellindex(k)+1));
    if prev == prev0
        ik = ik + 1;
    else
        ik = 1;
    end
    AO{1}.DeviceList(k,:) = [prev ik];
    AO{1}.ElementList(k,:)= str2num(map{k}(elemindex(k):sep(k)));
    prev0 = prev;
end

AO{1}.DeviceName(:,:)               = dev;
AO{1}.Monitor.TangoNames(:,:)       = strcat(dev, '/XPosDD');
AO{1}.CommonNames(:,:)              = [repmat('BPMx',nb,1) num2str((1:nb)','%02d')];

AO{1}.Status                        = ones(nb,1);
% AO{1}.Status([5:22]) = 0;
AO{1}.Monitor.HW2PhysicsParams(:,:) = 1e-3*ones(nb,1);
AO{1}.Monitor.Physics2HWParams(:,:) = 1e3*ones(nb,1);

% 2 lignes ajoutes pour test debug bpm versus le reste du monde
AO{1}.Monitor.Handles(:,1)       = NaN*ones(nb,1);
AO{1}.Monitor.DataType         = 'Vector';


if ControlRoomFlag
        AO{1}.GroupId = tango_group_create2('BPM');
        tango_group_add(AO{1}.GroupId,AO{1}.DeviceName');
else
   AO{1}.GroupId = NaN;
end


AO{1}.Status = AO{1}.Status(:);
AO{1}.Monitor.Mode            = Mode; %'Special';
AO{1}.Monitor.SpecialFunction = 'get_horbit_boo';

% BPM with 0 status
%AO{1}.Status([3 5 8 12]) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BPMz data: status field designates if BPM in use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vertical plane
AO{2} = AO{1};
AO{2}.FamilyName              = 'BPMz';
AO{2}.MemberOf                 = {'BPM'; 'VBPM'; 'Diagnostics'};
AO{2}.Monitor.TangoNames(:,:) = strcat(dev,'/ZPosDD');
AO{2}.CommonNames(:,:) = [repmat('BPMz',nb,1) num2str((1:nb)','%02d')];

% 2 lignes ajoutes pour test debug bpm versus le reste du monde
AO{2}.Monitor.DataType         = 'Vector';
AO{2}.Monitor.Handles(:,1)       = NaN*ones(nb,1);

AO{2}.Status = AO{2}.Status(:);

AO{2}.Monitor.Mode            = Mode; %'Special';
AO{2}.Monitor.SpecialFunction = 'get_vorbit_boo';

% BPM with 0 status
%AO{1}.Status([3 5 8 12]) = 0;

% %===========================================================
% %Corrector data: status field designates if corrector in use
% %===========================================================
% 
setao(cell2field(AO));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SLOW HORIZONTAL CORRECTORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ifam = 3;
AO{ifam}.FamilyName               = 'HCOR';
AO{ifam}.FamilyType               = 'COR';
AO{ifam}.MemberOf                 = {'MachineConfig'; 'COR'; 'MCOR'; 'HCOR'; 'Magnet'};

AO{ifam}.Monitor.Mode             = Mode;
AO{ifam}.Monitor.DataType         = 'Scalar';
AO{ifam}.Monitor.Units            = 'Hardware';
AO{ifam}.Monitor.HWUnits          = 'A';
AO{ifam}.Monitor.PhysicsUnits     = 'rad';
AO{ifam}.Monitor.HW2PhysicsFcn    = @amp2k;
AO{ifam}.Monitor.Physics2HWFcn    = @k2amp;


% Get mapping from TANGO static database
if ControlRoomFlag
    map = tango_get_db_property('booster','tracy_correctorH_mapping');
else
    map= {  'HCOR01::BOO-C01/AE/CH.01'
            'HCOR02::BOO-C02/AE/CH.02'
            'HCOR03::BOO-C03/AE/CH.03'
            'HCOR04::BOO-C04/AE/CH.04'
            'HCOR05::BOO-C05/AE/CH.05'
            'HCOR06::BOO-C06/AE/CH.06'
            'HCOR07::BOO-C07/AE/CH.07'
            'HCOR08::BOO-C08/AE/CH.08'
            'HCOR09::BOO-C09/AE/CH.09'
            'HCOR10::BOO-C10/AE/CH.10'
            'HCOR11::BOO-C11/AE/CH.11'
            'HCOR12::BOO-C12/AE/CH.12'
            'HCOR13::BOO-C13/AE/CH.13'
            'HCOR14::BOO-C14/AE/CH.14'
            'HCOR15::BOO-C15/AE/CH.15'
            'HCOR16::BOO-C16/AE/CH.16'
            'HCOR17::BOO-C17/AE/CH.17'
            'HCOR18::BOO-C18/AE/CH.18'
            'HCOR19::BOO-C19/AE/CH.19'
            'HCOR20::BOO-C20/AE/CH.20'
            'HCOR21::BOO-C21/AE/CH.21'
            'HCOR22::BOO-C22/AE/CH.22'};
    map=map';    
end    
dev = regexprep(map,'^HCOR\d*::','')';
cellindex = cell2mat(regexpi(map,'C[0-9]','once'))+1;
numindex  = cell2mat(regexpi(map,'\.[0-9]','once'))+1;
elemindex = cell2mat(regexpi(map,'HCOR[0-9]','once'))+4;
sep       = cell2mat(regexpi(map,'::','once'))-1;

nb = size(map,2);
for k = 1:nb
    AO{ifam}.DeviceList(k,:)  = [str2num(map{k}(cellindex(k):cellindex(k)+1)) 1];
    AO{ifam}.ElementList(k,:) = str2num(map{k}(elemindex(k):sep(k)));
end

AO{ifam}.DeviceName(:,:)            = dev;
AO{ifam}.CommonNames(:,:)           = [repmat(AO{ifam}.FamilyName,nb,1) num2str((1:nb)','%03d')];
AO{ifam}.Status                     = ones(nb,1);
AO{ifam}.Monitor.TangoNames(:,:)    = strcat(dev,'/current');
AO{ifam}.Monitor.Handles(:,1)       = NaN*ones(nb,1);

%Load fields from datablock
% AT use the "A-coefficients" for correctors plus an offset
[C, Leff, MagnetType, coefficients] = magnetcoefficients(AO{ifam}.FamilyName);

for ii=1:nb
    AO{ifam}.Monitor.HW2PhysicsParams{1}(ii,:)  = coefficients;
    AO{ifam}.Monitor.Physics2HWParams{1}(ii,:)  = coefficients;
end

AO{ifam}.Setpoint = AO{ifam}.Monitor;
AO{ifam}.Desired  = AO{ifam}.Monitor;

AO{ifam}.Setpoint.Range(:,:)        = repmat([-1.5 1.5],nb,1); % 1 A for ???
AO{ifam}.Setpoint.Tolerance(:,:)    = 1000*ones(nb,1);
AO{ifam}.Setpoint.DeltaRespMat(:,:) = ones(nb,1)*2e-5; % ???

AO{ifam}.Status = AO{ifam}.Status(:);
%convert response matrix kicks to HWUnits (after AO is loaded to AppData)
setao(cell2field(AO));   %required to make physics2hw function
AO{ifam}.Setpoint.DeltaRespMat = physics2hw(AO{ifam}.FamilyName,'Setpoint', ...
    AO{ifam}.Setpoint.DeltaRespMat, AO{ifam}.DeviceList);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SLOW VERTICAL CORRECTORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ifam = ifam +1;
AO{ifam}.FamilyName               = 'VCOR';
AO{ifam}.FamilyType               = 'COR';
AO{ifam}.MemberOf                 = {'MachineConfig'; 'COR'; AO{ifam}.FamilyName; 'Magnet'};

AO{ifam}.Monitor.Mode             = Mode;
AO{ifam}.Monitor.DataType         = 'Scalar';
AO{ifam}.Monitor.Units            = 'Hardware';
AO{ifam}.Monitor.HWUnits          = 'A';
AO{ifam}.Monitor.PhysicsUnits     = 'rad';
AO{ifam}.Monitor.HW2PhysicsFcn = @amp2k;
AO{ifam}.Monitor.Physics2HWFcn = @k2amp;

% Get mapping from TANGO static database
if ControlRoomFlag
    map = tango_get_db_property('booster','tracy_correctorV_mapping');
else
    map={   'VCOR01::BOO-C01/AE/CV.01'
            'VCOR02::BOO-C02/AE/CV.02'
            'VCOR03::BOO-C03/AE/CV.03'
            'VCOR04::BOO-C04/AE/CV.04'
            'VCOR05::BOO-C05/AE/CV.05'
            'VCOR06::BOO-C06/AE/CV.06'
            'VCOR07::BOO-C07/AE/CV.07'
            'VCOR08::BOO-C08/AE/CV.08'
            'VCOR09::BOO-C09/AE/CV.09'
            'VCOR10::BOO-C10/AE/CV.10'
            'VCOR11::BOO-C11/AE/CV.11'
            'VCOR12::BOO-C12/AE/CV.12'
            'VCOR13::BOO-C13/AE/CV.13'
            'VCOR14::BOO-C14/AE/CV.14'
            'VCOR15::BOO-C15/AE/CV.15'
            'VCOR16::BOO-C16/AE/CV.16'
            'VCOR17::BOO-C17/AE/CV.17'
            'VCOR18::BOO-C18/AE/CV.18'
            'VCOR19::BOO-C19/AE/CV.19'
            'VCOR20::BOO-C20/AE/CV.20'
            'VCOR21::BOO-C21/AE/CV.21'
            'VCOR22::BOO-C22/AE/CV.22'};
        map=map';
end    
dev = regexprep(map,'^VCOR\d*::','')';
cellindex = cell2mat(regexpi(map,'C[0-9]','once'))+1;
numindex  = cell2mat(regexpi(map,'\.[0-9]','once'))+1;
elemindex = cell2mat(regexpi(map,'VCOR[0-9]','once'))+4;
sep       = cell2mat(regexpi(map,'::','once'))-1;

nb = size(map,2);
for k = 1:nb
    AO{ifam}.DeviceList(k,:)  = [str2num(map{k}(cellindex(k):cellindex(k)+1)) 1];
    AO{ifam}.ElementList(k,:) = str2num(map{k}(elemindex(k):sep(k)));
end

AO{ifam}.DeviceName(:,:)            = dev;
AO{ifam}.Status                     = ones(nb,1);
AO{ifam}.CommonNames(:,:)           = [repmat(AO{ifam}.FamilyName,nb,1) num2str((1:nb)','%03d')];
AO{ifam}.Monitor.TangoNames(:,:)    = strcat(dev,'/current');

AO{ifam}.Monitor.Handles(:,1)    = NaN*ones(nb,1);

%Load fields from datablock
% AT use the "A-coefficients" for correctors plus an offset
[C, Leff, MagnetType, coefficients] = magnetcoefficients(AO{ifam}.FamilyName);

for ii = 1:nb
    AO{ifam}.Monitor.HW2PhysicsParams{1}(ii,:)  = coefficients;
    AO{ifam}.Monitor.Physics2HWParams{1}(ii,:)  = coefficients;
end

AO{ifam}.Setpoint = AO{ifam}.Monitor;
AO{ifam}.Desired  = AO{ifam}.Monitor;

AO{ifam}.Setpoint.Range(:,:) = repmat([-1.5 1.5],nb,1); %  1 A for ???
AO{ifam}.Setpoint.Tolerance(:,:) = 1000*ones(nb,1);
AO{ifam}.Setpoint.DeltaRespMat(:,:) = ones(nb,1)*4e-5; % ???

%convert response matrix kicks to HWUnits (after AO is loaded to AppData)
setao(cell2field(AO));   %required to make physics2hw function
AO{ifam}.Setpoint.DeltaRespMat = physics2hw(AO{ifam}.FamilyName,'Setpoint', ...
    AO{ifam}.Setpoint.DeltaRespMat, AO{ifam}.DeviceList);

%=============================
% MAIN MAGNETS
%=============================

%===========
%% Dipole data
%===========

% *** BEND ***
ifam = ifam+1;
AO{ifam}.FamilyName                 = 'BEND';
AO{ifam}.FamilyType                 = 'BEND';
AO{ifam}.MemberOf                   = {'MachineConfig'; 'BEND'; 'Magnet';};
HW2PhysicsParams                    = magnetcoefficients('BEND');
Physics2HWParams                    = HW2PhysicsParams;

AO{ifam}.Monitor.Mode               = Mode;
AO{ifam}.Monitor.DataType           = 'Scalar';
AO{ifam}.Monitor.Units              = 'Hardware';
AO{ifam}.Monitor.HW2PhysicsFcn      = @bend2gev;
AO{ifam}.Monitor.Physics2HWFcn      = @gev2bend;
AO{ifam}.Monitor.HWUnits            = 'A';
AO{ifam}.Monitor.PhysicsUnits       = 'energy';


% nb = size(map,2);
nb = 1;

AO{ifam}.DeviceList(:,:)   = [1 1];
AO{ifam}.ElementList(:,:)  = 1;
AO{ifam}.DeviceName(:,:)   = 'BOO/AE/Dipole';
AO{ifam}.Monitor.TangoNames(:,:)  = strcat(AO{ifam}.DeviceName,'/d1current');
AO{ifam}.OffsetR.TangoNames(:,:)  = strcat(AO{ifam}.DeviceName,'/d1waveformOffset');

AO{ifam}.Status = 1;
AO{ifam}.Monitor.Handles = NaN;

val = 1.0;
AO{ifam}.Monitor.HW2PhysicsParams{1}(1,:)                 = HW2PhysicsParams;
AO{ifam}.Monitor.HW2PhysicsParams{2}(1,:)                 = val;
AO{ifam}.Monitor.Physics2HWParams{1}(1,:)                 = Physics2HWParams;
AO{ifam}.Monitor.Physics2HWParams{2}(1,:)                 = val;

% When using bend2gev / gev2bend don't use Params
AO{ifam}.Monitor = rmfield(AO{ifam}.Monitor,'HW2PhysicsParams');
AO{ifam}.Monitor = rmfield(AO{ifam}.Monitor,'Physics2HWParams');

AO{ifam}.Setpoint = AO{ifam}.Monitor;
AO{ifam}.Desired  = AO{ifam}.Monitor;
AO{ifam}.OffsetR  = AO{ifam}.Monitor;
AO{ifam}.OffsetW  = AO{ifam}.Monitor;
AO{ifam}.Setpoint.TangoNames = strcat(AO{ifam}.DeviceName,'/current');
AO{ifam}.OffsetW.TangoNames = strcat(AO{ifam}.DeviceName,'/waveformOffset');

AO{ifam}.Setpoint.Range(:,:) = [0 560]; % 525 A for 1.71T
AO{ifam}.Setpoint.Tolerance(:,:) = 0.05;
AO{ifam}.Setpoint.DeltaRespMat(:,:) = 0.05;

%convert response matrix kicks to HWUnits (after AO is loaded to AppData)
%setao(cell2field(AO));   %required to make physics2hw function
%AO{ifam}.Setpoint.DeltaRespMat=physics2hw(AO{ifam}.FamilyName,'Setpoint',AO{ifam}.Setpoint.DeltaRespMat,AO{ifam}.DeviceList);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QUADRUPOLE MAGNETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get mapping from TANGO static database
if ControlRoomFlag
    map = tango_get_db_property('booster','tracy_quadrupole_mapping');
else
    map={   'QD::BOO/AE/QD'
            'QF::BOO/AE/QF'};
    map=map';    
end    
iifam   = cell2mat(regexpi(map,'/Q','once'))+2;
ifound  = regexpi(map,'/QF');

%%% build mapping for the two families
%%% cleanup
mapQF = []; mapQD = [];
for k = 1:length(ifound)
    if isempty(ifound{k})
        mapQD = [mapQD map(k)];
    else
        mapQF = [mapQF map(k)];
    end
end


mapQP  = {mapQF,mapQD};
nameQP = {'QF', 'QD'};

% for both family
for k = 1:2    
    sep = cell2mat(regexprep(mapQP{k},'::','once'))-1;
    dev = regexprep(mapQP{k},'^Q[F,D]\d*::','')';
    ifam = ifam + 1;

    AO{ifam}.FamilyName                 = nameQP{k};
    AO{ifam}.FamilyType                 = 'QUAD';
    AO{ifam}.MemberOf                   = {'MachineConfig'; 'QUAD'; 'Magnet';};
    HW2PhysicsParams                    = magnetcoefficients(AO{ifam}.FamilyName);
    Physics2HWParams                    = magnetcoefficients(AO{ifam}.FamilyName);

    AO{ifam}.Monitor.Mode               = Mode;
    AO{ifam}.Monitor.DataType           = 'Scalar';
    AO{ifam}.Monitor.Units              = 'Hardware';
    AO{ifam}.Monitor.HWUnits            = 'A';
    AO{ifam}.Monitor.PhysicsUnits       = 'meter^-2';
    AO{ifam}.Monitor.HW2PhysicsFcn      = @amp2k;
    AO{ifam}.Monitor.Physics2HWFcn      = @k2amp;
    
    AO{ifam}.DeviceList(1,:) = [1 1];
    AO{ifam}.ElementList(1,:)= 1;

    prev0 = prev;

    AO{ifam}.DeviceName(:,:)    = dev;
    
    %% Build common names
    AO{ifam}.CommonNames(:,:) = AO{ifam}.FamilyName;

    AO{ifam}.Monitor.TangoNames(:,:)  = strcat(dev,'/current');

    AO{ifam}.Status = 1;
    AO{ifam}.Monitor.Handles(:,1) = NaN;

    AO{ifam}.Monitor.HW2PhysicsParams{1}(1,:)                 = HW2PhysicsParams;
    AO{ifam}.Monitor.HW2PhysicsParams{2}(1,:)                 = val;
    AO{ifam}.Monitor.Physics2HWParams{1}(1,:)                 = Physics2HWParams;
    AO{ifam}.Monitor.Physics2HWParams{2}(1,:)                 = val;

    AO{ifam}.Setpoint = AO{ifam}.Monitor;
    AO{ifam}.Desired  = AO{ifam}.Monitor;

    AO{ifam}.Setpoint.Range(:,:) = [0 300]; %260 A for 19.7Tm-1
    AO{ifam}.Setpoint.Tolerance(:,:) = 0.05;
    AO{ifam}.Setpoint.DeltaRespMat(:,:) = 0.02;

    %convert response matrix kicks to HWUnits (after AO is loaded to AppData)
    setao(cell2field(AO));   %required to make physics2hw function
%     AO{ifam}.Setpoint.DeltaRespMat = physics2hw(AO{ifam}.FamilyName, ...
%         'Setpoint',AO{ifam}.Setpoint.DeltaRespMat,AO{ifam}.DeviceList);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SEXTUPOLE MAGNETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get mapping from TANGO static database
if ControlRoomFlag
    map=tango_get_db_property('booster','tracy_sextupole_mapping');
else
    map={   'SD::BOO/AE/SD'
            'SF::BOO/AE/SF'};        
    map=map';
end
iifam = cell2mat(regexpi(map,'/S','once'))+2;
ifound  = regexpi(map,'/SF');

%%% build mapping for the ten families
%%% cleanup
mapSF = []; mapSD = [];
for k = 1:length(ifound)
    if isempty(ifound{k})
        mapSD = [mapSD map(k)];
    else
        mapSF = [mapSF map(k)];
    end
end


mapS  = {mapSF,mapSD};
nameS = {'SF', 'SD'};

% For both families
for k = 1:2
    sep = cell2mat(regexpi(mapS{k},'::','once'))-1;
    dev = regexprep(mapS{k},'^S[F,D]\d*::','')';

    ifam=ifam+1;

    AO{ifam}.FamilyName                = nameS{k};
    AO{ifam}.FamilyType                = 'SEXT';
    AO{ifam}.MemberOf                  = {'MachineConfig'; 'SEXT'; 'Magnet';};
    HW2PhysicsParams                   = magnetcoefficients(AO{ifam}.FamilyName);
    Physics2HWParams                   = magnetcoefficients(AO{ifam}.FamilyName);

    AO{ifam}.Monitor.Mode              = Mode;
    AO{ifam}.Monitor.DataType          = 'Scalar';
    AO{ifam}.Monitor.Units             = 'Hardware';
    AO{ifam}.Monitor.HW2PhysicsFcn     = @amp2k;
    AO{ifam}.Monitor.Physics2HWFcn     = @k2amp;
    AO{ifam}.Monitor.HWUnits           = 'A';
    AO{ifam}.Monitor.PhysicsUnits      = 'meter^-3';

    AO{ifam}.DeviceList(:,:) = [1 1];
    AO{ifam}.ElementList(:,:)= 1;

    AO{ifam}.DeviceName(:,:)    = dev;
    %% Build common names
    AO{ifam}.CommonNames(:,:)   = AO{ifam}.FamilyName;
    AO{ifam}.Monitor.TangoNames(:,:)  = strcat(dev,'/current');
    

    AO{ifam}.Status = 1;
    AO{ifam}.Monitor.Handles(:,1) = NaN;

    val = 1.0;
    AO{ifam}.Monitor.HW2PhysicsParams{1}(1,:)               = HW2PhysicsParams;
    AO{ifam}.Monitor.HW2PhysicsParams{2}(1,:)               = val;
    AO{ifam}.Monitor.Physics2HWParams{1}(1,:)               = Physics2HWParams;
    AO{ifam}.Monitor.Physics2HWParams{2}(1,:)               = val;

    AO{ifam}.Setpoint = AO{ifam}.Monitor;
    AO{ifam}.Desired  = AO{ifam}.Monitor;

    AO{ifam}.Setpoint.Range(:,:) = [0 400]; % 343 A for 320 Tm-2
    AO{ifam}.Setpoint.Tolerance(:,:) = 0.05;
    AO{ifam}.Setpoint.DeltaRespMat(:,:) = 0.05;

    %convert response matrix kicks to HWUnits (after AO is loaded to AppData)
    setao(cell2field(AO));   %required to make physics2hw function
%     AO{ifam}.Setpoint.DeltaRespMat=physics2hw(AO{ifam}.FamilyName,'Setpoint',AO{ifam}.Setpoint.DeltaRespMat,AO{ifam}.DeviceList);
end

% Convert to new format
AO = cell2field(AO);


%=======
%% TUNE
%=======
ifam = 'TUNE';
AO.(ifam).FamilyName  = ifam;
AO.(ifam).FamilyType  = 'Diagnostic';
AO.(ifam).MemberOf    = {'Diagnostics'};
AO.(ifam).CommonNames = ['nux';'nuz';'nus'];
AO.(ifam).DeviceList  = [1 1; 1 2; 1 3];
AO.(ifam).ElementList = [1 2 3]';
AO.(ifam).Status      = [1 1 1]';

AO.(ifam).Monitor.Mode                   = Mode; %'Simulator';  % Mode;
AO.(ifam).Monitor.DataType               = 'Scalar';
AO.(ifam).Monitor.DataTypeIndex          = [1 2];
AO.(ifam).Monitor.TangoNames              = ['BOO/DG/BPM-TUNEX/Nu'; ...
                                             'BOO/DG/BPM-TUNEZ/Nu'; 'BOO/DG/BPM-TUNEX/Nu'];
AO.(ifam).Monitor.Units                  = 'Hardware';
AO.(ifam).Monitor.Handles                = NaN;
AO.(ifam).Monitor.HW2PhysicsParams       = 1;
AO.(ifam).Monitor.Physics2HWParams       = 1;
AO.(ifam).Monitor.HWUnits                = 'fractional tune';
AO.(ifam).Monitor.PhysicsUnits           = 'fractional tune';

%%%%%%%%%%%%%%%%%%
%% Pulsed Magnet
%%%%%%%%%%%%%%%%%%

%% Injection kicker
ifam = 'K_Inj'; 
AO.(ifam).FamilyName           = 'K_Inj';
AO.(ifam).MemberOf             = {'Injection';'Archivable';'EP'};
AO.(ifam).Monitor.Mode         = Mode;
AO.(ifam).Monitor.DataType     = 'Scalar';

AO.(ifam).Status                   = 1;
AO.(ifam).DeviceName               = cellstr('BOO-C01/EP/AL_K.Inj');
AO.(ifam).CommonNames              = 'K_Inj';
AO.(ifam).ElementList              = 1;
AO.(ifam).DeviceList(:,:)          = [1 1];
AO.(ifam).Monitor.Handles(:,1)     = NaN;
AO.(ifam).Monitor.TangoNames       = strcat(AO.(ifam).DeviceName, '/voltage');
AO.(ifam).Monitor.HW2PhysicsParams = 1;
AO.(ifam).Monitor.Physics2HWParams = 1;
AO.(ifam).Monitor.Units            = 'Hardware';   
AO.(ifam).Monitor.HWUnits          = 'V';
AO.(ifam).Monitor.PhysicsUnits     = 'mrad';
AO.(ifam).Setpoint = AO.(ifam).Monitor;
AO.(ifam).Desired = AO.(ifam).Monitor;

%% Injection Septum
ifam = 'SEP_Inj'; 
AO.(ifam).FamilyName           = 'SEP_Inj';
AO.(ifam).MemberOf             = {'Injection';'Archivable';'EP'};
AO.(ifam).Monitor.Mode         = Mode;
AO.(ifam).Monitor.DataType     = 'Scalar';

AO.(ifam).Status                   = 1;
AO.(ifam).DeviceName               = cellstr('BOO-C22/EP/AL_SEP_P.Inj');
AO.(ifam).CommonNames              = 'SEP_P';
AO.(ifam).ElementList              = 1;
AO.(ifam).DeviceList(:,:)           = [1 1];
AO.(ifam).Monitor.Handles(:,1)     = NaN;
AO.(ifam).Monitor.TangoNames       = strcat(AO.(ifam).DeviceName, '/voltage');
AO.(ifam).Monitor.HW2PhysicsParams = 1;
AO.(ifam).Monitor.Physics2HWParams = 1;
AO.(ifam).Monitor.Units            = 'Hardware';   
AO.(ifam).Monitor.HWUnits          = 'V';
AO.(ifam).Monitor.PhysicsUnits     = 'mrad';
AO.(ifam).Setpoint = AO.(ifam).Monitor;
AO.(ifam).Desired = AO.(ifam).Monitor;

%% Injection Septum
ifam = 'SEPA_EXT'; 
AO.(ifam).FamilyName           = ifam;
AO.(ifam).MemberOf             = {'Injection';'Archivable';'EP'};
AO.(ifam).Monitor.Mode         = Mode;
AO.(ifam).Monitor.DataType     = 'Scalar';

AO.(ifam).Status                   = 0;
AO.(ifam).DeviceName               = cellstr('BOO-C22/EP/AL_SEP_P.Inj');
AO.(ifam).CommonNames              = 'SEP_P';
AO.(ifam).ElementList              = 1;
AO.(ifam).DeviceList(:,:)           = [1 1];
AO.(ifam).Monitor.Handles(:,1)     = NaN;
AO.(ifam).Monitor.TangoNames       = strcat(AO.(ifam).DeviceName, '/voltage');
AO.(ifam).Monitor.HW2PhysicsParams = 1;
AO.(ifam).Monitor.Physics2HWParams = 1;
AO.(ifam).Monitor.Units            = 'Hardware';   
AO.(ifam).Monitor.HWUnits          = 'V';
AO.(ifam).Monitor.PhysicsUnits     = 'mrad';
AO.(ifam).Setpoint = AO.(ifam).Monitor;
AO.(ifam).Desired = AO.(ifam).Monitor;

%% Injection Septum
ifam = 'SEP_P_EXT'; 
AO.(ifam).FamilyName           = ifam;
AO.(ifam).MemberOf             = {'Injection';'Archivable';'EP'};
AO.(ifam).Monitor.Mode         = Mode;
AO.(ifam).Monitor.DataType     = 'Scalar';

AO.(ifam).Status                   = 0;
AO.(ifam).DeviceName               = cellstr('BOO-C11/EP/AL_SEP_P.Ext');
AO.(ifam).CommonNames              = 'SEP_P_EXT';
AO.(ifam).ElementList              = 1;
AO.(ifam).DeviceList(:,:)           = [1 1];
AO.(ifam).Monitor.Handles(:,1)     = NaN;
AO.(ifam).Monitor.TangoNames       = strcat(AO.(ifam).DeviceName, '/voltage');
AO.(ifam).Monitor.HW2PhysicsParams = 1;
AO.(ifam).Monitor.Physics2HWParams = 1;
AO.(ifam).Monitor.Units            = 'Hardware';   
AO.(ifam).Monitor.HWUnits          = 'V';
AO.(ifam).Monitor.PhysicsUnits     = 'mrad';
AO.(ifam).Setpoint = AO.(ifam).Monitor;
AO.(ifam).Desired = AO.(ifam).Monitor;

%%%%%%%%%%%%%%%%%%
%% VACUUM SYSTEM
%%%%%%%%%%%%%%%%%%

%% IonPump
ifam = 'PI';
AO.(ifam).FamilyName           = 'PI';
AO.(ifam).FamilyType           = 'PI';
AO.(ifam).MemberOf             = {'PlotFamily'; 'IonPump'; 'Pressure'; 'Archivable'};
AO.(ifam).Monitor.Mode         = Mode;
AO.(ifam).Monitor.DataType     = 'Scalar';
if ControlRoomFlag
    map     = tango_get_db_property('booster','pompe_ionique');
else
    map={'BOO-B1/VI/PI55.1'
        'BOO-B1/VI/PI55.2'
        'BOO-B1/VI/PI55.3'
        'BOO-B1/VI/PI55.4'
        'BOO-B1/VI/PI55.5'
        'BOO-B1/VI/PI55.6'
        'BOO-B1/VI/PI55.7'
        'BOO-B1/VI/PI55.8'
        'BOO-B1/VI/PI55.9'
        'BOO-Binj/VI/PI55.1'
        'BOO-Binj/VI/PI55.2'
        'BOO-Binj/VI/PI300.3'
        'BOO-Binj/VI/PI55.4'
        'BOO-Binj/VI/PI55.5'
        'BOO-Binj/VI/PI55.6'
        'BOO-Binj/VI/PI500.7'
        'BOO-B2/VI/PI55.1'
        'BOO-B2/VI/PI55.2'
        'BOO-B2/VI/PI55.3'
        'BOO-B2/VI/PI55.4'
        'BOO-B2/VI/PI55.5'
        'BOO-B2/VI/PI55.6'
        'BOO-B2/VI/PI55.7'
        'BOO-Bext/VI/PI55.1'
        'BOO-Bext/VI/PI55.2'
        'BOO-Bext/VI/PI300.3'
        'BOO-Bext/VI/PI55.4'
        'BOO-Bext/VI/PI55.5'
        'BOO-Bext/VI/PI55.6'
        'BOO-B3/VI/PI55.1'
        'BOO-B3/VI/PI55.2'
        'BOO-B3/VI/PI55.3'};
map=map';
end    
numindex  = cell2mat(regexpi(map,'\.[0-9]','once'))+1;

nb = size(map,2);
for k = 1:nb,
    switch map{k}(1:6)
        case 'BOO-B1'
            k1 = 1;
        case 'BOO-Bi'
            k1 = 2;
        case 'BOO-B2'
            k1 = 3;
        case 'BOO-Be'
            k1 = 4;
        case 'BOO-B3'
            k1 = 5;
        case 'BOO-B4'
            k1 = 6;
    end            
    AO.(ifam).DeviceList(k,:)  = [k1 str2double(map{k}(numindex(k)))];
    AO.(ifam).ElementList(k,:) = k;
    AO.(ifam).CommonNames(k,:) = ['PI' num2str(ik)];
end

AO.(ifam).Status                   = ones(nb,1);
%AO.(ifam).Status(26) = 0;% BOO-Bext/VI/PI300.3
AO.(ifam).DeviceName               = map';
AO.(ifam).CommonNames              = AO.(ifam).CommonNames';
AO.(ifam).ElementList              = (1:nb)';
AO.(ifam).Monitor.Handles(:,1)     = NaN*ones(nb,1);
AO.(ifam).Monitor.TangoNames       = strcat(AO.(ifam).DeviceName, '/pressure');
AO.(ifam).Monitor.HW2PhysicsParams = 1;
AO.(ifam).Monitor.Physics2HWParams = 1;
AO.(ifam).Monitor.Units            = 'Hardware';   
AO.(ifam).Monitor.HWUnits          = 'mBar';
AO.(ifam).Monitor.PhysicsUnits     = 'mBar';

%% Penning Gauges
ifam = 'JPEN';
AO.(ifam).FamilyName           = 'JPEN';
AO.(ifam).MemberOf             = {'PlotFamily'; 'PenningGauge'; 'Pressure'; 'Archivable'};
AO.(ifam).Monitor.Mode         = Mode;
AO.(ifam).Monitor.DataType     = 'Scalar';

if ControlRoomFlag
    map     = tango_get_db_property('booster','jauge_penning');
else
    map={'BOO-B1/VI/JPEN.1'
        'BOO-B1/VI/JPEN.2'
        'BOO-Binj/VI/JPEN.1'
        'BOO-Binj/VI/JPEN.2'
        'BOO-B2/VI/JPEN.1'
        'BOO-B2/VI/JPEN.2'
        'BOO-B2/VI/JPEN.3'
        'BOO-Bext/VI/JPEN.1'
        'BOO-Bext/VI/JPEN.2'
        'BOO-B3/VI/JPEN.1'
        'BOO-B3/VI/JPEN.2'
        'BOO-B4/VI/JPEN.1'
        'BOO-B4/VI/JPEN.2'
        'BOO-B4/VI/JPEN.3'};
        map=map';
end
numindex  = cell2mat(regexpi(map,'\.[0-9]','once'))+1;

nb = size(map,2);
for k = 1:nb,
    switch map{k}(1:6)
        case 'BOO-B1'
            k1 = 1;
        case 'BOO-Bi'
            k1 = 2;
        case 'BOO-B2'
            k1 = 3;
        case 'BOO-Be'
            k1 = 4;
        case 'BOO-B3'
            k1 = 5;
        case 'BOO-B4'
            k1 = 6;
    end            
    AO.(ifam).DeviceList(k,:)  = [k1 str2double(map{k}(numindex(k)))];
    AO.(ifam).ElementList(k,:) = k;
    AO.(ifam).CommonNames(k,:) = ['JPEN' num2str(ik)];
end

AO.(ifam).Status                   = ones(nb,1);
AO.(ifam).DeviceName               = map';
AO.(ifam).CommonNames              = AO.(ifam).CommonNames';
AO.(ifam).ElementList              = (1:nb)';
AO.(ifam).Monitor.Handles(:,1)     = NaN*ones(nb,1);
AO.(ifam).Monitor.TangoNames       = strcat(AO.(ifam).DeviceName, '/pressure');
AO.(ifam).Monitor.HW2PhysicsParams = 1;
AO.(ifam).Monitor.Physics2HWParams = 1;
AO.(ifam).Monitor.Units            = 'Hardware';   
AO.(ifam).Monitor.HWUnits          = 'mBar';
AO.(ifam).Monitor.PhysicsUnits     = 'mBar';

%% Pirani Gauges
ifam = 'JPIR';
AO.(ifam).FamilyName           = 'JPIR';
AO.(ifam).MemberOf             = {'PlotFamily'; 'PiraniGauge'; 'Pressure'; 'Archivable'};
AO.(ifam).Monitor.Mode         = Mode;
AO.(ifam).Monitor.DataType     = 'Scalar';

if ControlRoomFlag
    map     = tango_get_db_property('booster','jauge_pirani');
else
    map={'BOO-B1/VI/JPIR.1'
        'BOO-Binj/VI/JPIR.1'
        'BOO-Binj/VI/JPIR.2'
        'BOO-B2/VI/JPIR.1'
        'BOO-Bext/VI/JPIR.1'
        'BOO-B3/VI/JPIR.1'
        'BOO-B4/VI/JPIR.1'};
    map=map';
end

numindex  = cell2mat(regexpi(map,'\.[0-9]','once'))+1;

nb = size(map,2);
for k = 1:nb,
    switch map{k}(1:6)
        case 'BOO-B1'
            k1 = 1;
        case 'BOO-Bi'
            k1 = 2;
        case 'BOO-B2'
            k1 = 3;
        case 'BOO-Be'
            k1 = 4;
        case 'BOO-B3'
            k1 = 5;
        case 'BOO-B4'
            k1 = 6;
    end            
    AO.(ifam).DeviceList(k,:)  = [k1 str2double(map{k}(numindex(k)))];
    AO.(ifam).ElementList(k,:) = k;
    AO.(ifam).CommonNames(k,:) = ['JPIR' num2str(ik)];
end

AO.(ifam).Status                   = ones(nb,1);
AO.(ifam).DeviceName               = map';
AO.(ifam).CommonNames              = AO.(ifam).CommonNames';
AO.(ifam).ElementList              = (1:nb)';
AO.(ifam).Monitor.Handles(:,1)     = NaN*ones(nb,1);
AO.(ifam).Monitor.TangoNames       = strcat(AO.(ifam).DeviceName, '/pressure');
AO.(ifam).Monitor.HW2PhysicsParams = 1;
AO.(ifam).Monitor.Physics2HWParams = 1;
AO.(ifam).Monitor.Units            = 'Hardware';   
AO.(ifam).Monitor.HWUnits          = 'mBar';
AO.(ifam).Monitor.PhysicsUnits     = 'mBar';

%====
%% DCCT
%====
ifam=ifam+1;
AO.DCCT.FamilyName                     = 'DCCT';
AO.DCCT.FamilyType                     = 'Diagnostic';
AO.DCCT.MemberOf                       = {'Diagnostics','Archivable','Plotfamily'};
AO.DCCT.CommonNames                    = 'DCCT';
AO.DCCT.DeviceList                     = [1 1; 1 2];
AO.DCCT.ElementList                    = (1:2)';
AO.DCCT.Status                         = AO.DCCT.ElementList;

AO.DCCT.Monitor.Mode                   = Mode;
AO.DCCT.FamilyName                     = 'DCCT';
AO.DCCT.deviceName                     = 'BOO-C01/DG/DCCT';
AO.DCCT.Monitor.DataType               = 'Scalar';
AO.DCCT.Monitor.TangoNames             = ['BOO-C01/DG/DCCT/iInj'; 'BOO-C01/DG/DCCT/iExt']; %afin de ne pas avoir de bug
AO.DCCT.Monitor.Units                  = 'Hardware';
AO.DCCT.Monitor.Handles                = NaN;
AO.DCCT.Monitor.HWUnits                = 'milli-ampere';
AO.DCCT.Monitor.PhysicsUnits           = 'A';
AO.DCCT.Monitor.HW2PhysicsParams       = 1;
AO.DCCT.Monitor.Physics2HWParams       = 1;

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
AO.(ifam).DeviceName(:,:)           = {'BOO/RF/MasterClock/'};

%Frequency Readback
AO.(ifam).Monitor.Mode                = Mode;
AO.(ifam).Monitor.DataType            = 'Scalar';
AO.(ifam).Monitor.Units               = 'Hardware';
AO.(ifam).Monitor.HW2PhysicsParams    = 1e+6;       %no hw2physics function necessary
AO.(ifam).Monitor.Physics2HWParams    = 1e-6;
AO.(ifam).Monitor.HWUnits             = 'MHz';
AO.(ifam).Monitor.PhysicsUnits        = 'Hz';
AO.(ifam).Monitor.TangoNames          = {'ANS/RF/MasterClock/frequency'};
AO.(ifam).Monitor.Handles             = NaN;
AO.(ifam).Monitor.Range               = [350e6 360e6];

AO.(ifam).Setpoint = AO.(ifam).Monitor;
AO.(ifam).Desired  = AO.(ifam).Monitor;


%=====================
%% Gamma Monitors DOSE
%====================
ifam = ifam+1;
AO.CIGdose.FamilyName                     = 'CIGdose';
AO.CIGdose.FamilyType                     = 'Radioprotection';
AO.CIGdose.MemberOf                       = {'Radioprotection','Archivable','Plotfamily'};
AO.CIGdose.CommonNames                    = 'CIG';

if ControlRoomFlag
    map = tango_get_db_property('anneau','gammamonitor_mapping');
else
    map={'LIN/RP/CIG.001'
        'LT1/RP/CIG.002'
        'BOO-C21/RP/CIG.003'
        'BOO-C04/RP/CIG.004'
        'LT2/RP/CIG.006'
        'BOO-C16/RP/CIG.007'
        'BOO-C17/RP/CIG.008'
        'BOO-C13/RP/CIG.009'
        'ANS-C16/RP/CIG.010'
        'ANS-C16/RP/CIG.011'
        'ANS-C01/RP/CIG.012'
        'ANS-C02/RP/CIG.013'
        'ANS-C02/RP/CIG.014'
        'ANS-C03/RP/CIG.015'
        'ANS-C04/RP/CIG.016'
        'ANS-C05/RP/CIG.017'
        'ANS-C05/RP/CIG.018'
        'ANS-C06/RP/CIG.019'
        'ANS-C07/RP/CIG.020'
        'ANS-C08/RP/CIG.021'
        'ANS-C09/RP/CIG.022'
        'ANS-C10/RP/CIG.023'
        'ANS-C11/RP/CIG.024'
        'ANS-C12/RP/CIG.025'
        'ANS-C12/RP/CIG.026'
        'ANS-C13/RP/CIG.027'
        'ANS-C14/RP/CIG.028'
        'ANS-C15/RP/CIG.029'
        'ANS-C16/RP/CIG.030'
        'ANS-C01/RP/CIG.031'
        'ANS-C03/RP/CIG.032'
        'ANS-C04/RP/CIG.033'
        'ANS-C02/RP/CIG.034'
        'ANS-C03/RP/CIG.035'
        'BOO-C06/RP/CIG.036'
        'BOO-C06/RP/CIG.037'
        'BOO-C02/RP/CIG.038'
        'BOO-C03/RP/CIG.039'
        'D01-1/RP/CIG'
        'D03-1/RP/CIG'
        'I04-M/RP/CIG'
        'D05-1/RP/CIG'
        'I05-L/RP/CIG'
        'I06-C/RP/CIG'
        'I08-M/RP/CIG'
        'D09-1/RP/CIG'
        'I10-C/RP/CIG'
        'I11-C/RP/CIG'
        'D13-1/RP/CIG'
        'I15-M/RP/CIG'
        'I16-M/RP/CIG'};
    map=map';
end
AO.CIGdose.DeviceName = map';

AO.CIGdose.Monitor.Mode                   = Mode;
AO.CIGdose.FamilyName                     = 'CIGdose';
AO.CIGdose.Monitor.DataType               = 'Scalar';

nb = length(AO.CIGdose.DeviceName); 
AO.CIGdose.DeviceList                     = [ones(1,nb); (1:nb)]';
AO.CIGdose.ElementList                    = (1:nb)';
AO.CIGdose.Status                         = ones(nb,1);

AO.CIGdose.Monitor.TangoNames            = strcat(AO.CIGdose.DeviceName,'/dose'); 

;%afin de ne pas avoir de bug
AO.CIGdose.Monitor.Units                  = 'Hardware';
AO.CIGdose.Monitor.Handles                = NaN;
AO.CIGdose.Monitor.HWUnits                = 'uGy';
AO.CIGdose.Monitor.PhysicsUnits           = 'uGy';
AO.CIGdose.Monitor.HW2PhysicsParams       = 1;
AO.CIGdose.Monitor.Physics2HWParams       = 1;

%=========================
%% Gamma Monitors DOSErate
%=========================

ifam = ifam+1;

AO.CIGrate = AO.CIGdose;
AO.CIGrate.FamilyName                     = 'CIGrate';
AO.CIGrate.CommonNames                    = 'CIGrate';
AO.CIGrate.FamilyName                     = 'CIGrate';

nb = length(AO.CIGrate.DeviceName); 
AO.CIGrate.Status                         = ones(nb,1);

AO.CIGrate.Monitor.TangoNames            = strcat(AO.CIGrate.DeviceName,'/doseRate'); 

%afin de ne pas avoir de bug
AO.CIG.Monitor.HWUnits                = 'uGy/h';
AO.CIG.Monitor.PhysicsUnits           = 'uGy/h';

%==================
%Machine Parameters
%==================
AO.MachineParameters.FamilyName                = 'MachineParameters';
AO.MachineParameters.FamilyType                = 'Parameter';
AO.MachineParameters.MemberOf                  = {'Diagnostics'};
AO.MachineParameters.Status                    = [1 1 1 1]';

AO.MachineParameters.Monitor.Mode              = Mode;
AO.MachineParameters.Monitor.DataType          = 'Scalar';
AO.MachineParameters.Monitor.Units             = 'Hardware';

%use spear2 process variable names
mp={
    'mode    '    'SPEAR:BeamStatus  '          [1 1]  1; ...
    'energy  '    'SPEAR:Energy      '          [1 2]  2; ...
    'current '    'SPEAR:BeamCurrAvg '          [1 3]  3; ...
    'lifetime'    'SPEAR:BeamLifetime'          [1 4]  4; ...
    };
AO.MachineParameters.Monitor.HWUnits          = ' ';
AO.MachineParameters.Monitor.PhysicsUnits     = ' ';

AO.MachineParameters.Setpoint.HWUnits         = ' ';
AO.MachineParameters.Setpoint.PhysicsUnits    = ' ';

for ii=1:size(mp,1)
    name  =mp(ii,1);    AO.MachineParameters.CommonNames(ii,:)            = char(name{1});
    %     name  =mp(ii,2);    AO.MachineParameters.Monitor.ChannelNames(ii,:)   = char(name{1});
    %     name  =mp(ii,2);    AO.MachineParameters.Setpoint.ChannelNames(ii,:)  = char(name{1});
    val   =mp(ii,3);    AO.MachineParameters.DeviceList(ii,:)             = val{1};
    val   =mp(ii,4);    AO.MachineParameters.ElementList(ii,:)            = val{1};

    AO.MachineParameters.Monitor.HW2PhysicsParams(ii,:)    = 1;
    AO.MachineParameters.Monitor.Physics2HWParams(ii,:)    = 1;
    AO.MachineParameters.Monitor.Handles(ii,1)  = NaN;
    AO.MachineParameters.Setpoint.HW2PhysicsParams(ii,:)   = 1;
    AO.MachineParameters.Setpoint.Physics2HWParams(ii,:)   = 1;
    AO.MachineParameters.Setpoint.Handles(ii,1)  = NaN;
end

% Save AO
setao(AO);

% The operational mode sets the path, filenames, and other important params
% Run setoperationalmode after most of the AO is built so that the Units and Mode fields
% can be set in setoperationalmode

setoperationalmode(OperationalMode);

%======================================================================
%======================================================================
%% Append Accelerator Toolbox information
%======================================================================
%======================================================================
disp('** Initializing Accelerator Toolbox information');

AT = getao;

ATindx = atindex(THERING);  %structure with fields containing indices

s = findspos(THERING,1:length(THERING)+1)';

% % Horizontal BPMS
% % WARNING: BPM1 is the one before the injection straigth section
% %          since a cell begins from begin of Straigths
% % CELL1 BPM1 to BPM7
% ifam = ('BPMx');
% AO.(ifam).AT.ATType  = ifam;
% AO.(ifam).AT.ATIndex = ATindx.BPM(:);
% AO.(ifam).Position   = s(AO.(ifam).AT.ATIndex);
% 
% % Vertical BPMS
% ifam = ('BPMz');
% AO.(ifam).AT.ATType  = ifam;
% AO.(ifam).AT.ATIndex = ATindx.BPM(:);
% AO.(ifam).Position   = s(AO.(ifam).AT.ATIndex);  
% 
% % SLOW HORIZONTAL CORRECTORS
% ifam = ('HCOR');
% AO.(ifam).AT.ATType  = ifam;
% AO.(ifam).AT.ATIndex = ATindx.HCOR(:);
% AO.(ifam).Position   = s(AO.(ifam).AT.ATIndex);
% 
% % SLOW VERTICAL CORRECTORS
% ifam = ('VCOR');
% AO.(ifam).AT.ATType  = ifam;
% AO.(ifam).AT.ATIndex = ATindx.VCOR(:);
% AO.(ifam).Position   = s(AO.(ifam).AT.ATIndex);  
% 
% %% BENDING magnets
% ifam = ('BEND');
% AO.(ifam).AT.ATType  = ifam;
% AO.(ifam).AT.ATIndex = ATindx.BEND(:);
% AO.(ifam).Position   = s(AO.(ifam).AT.ATIndex);
% 
% %% QF magnets
% ifam = ('QF');
% AO.(ifam).AT.ATType  = 'QUAD';
% AO.(ifam).AT.ATIndex = ATindx.QPF(:);
% AO.(ifam).Position   = s(AO.(ifam).AT.ATIndex);
% AO.(ifam).AT.ATParameterGroup = mkparamgroup(THERING,AO.(ifam).AT.ATIndex,'K');
% 
% %% QD magnets
% ifam = ('QD');
% AO.(ifam).AT.ATType  = 'QUAD';
% AO.(ifam).AT.ATIndex = ATindx.QPD(:);
% AO.(ifam).Position   = s(AO.(ifam).AT.ATIndex);
% AO.(ifam).AT.ATParameterGroup = mkparamgroup(THERING,AO.(ifam).AT.ATIndex,'K');
% 
% %% SF magnets
% ifam = ('SF');
% AO.(ifam).AT.ATType  = 'SEXT';
% AO.(ifam).AT.ATIndex = ATindx.SXF(:);
% AO.(ifam).Position   = s(AO.(ifam).AT.ATIndex);
% AO.(ifam).AT.ATParameterGroup = mkparamgroup(THERING,AO.(ifam).AT.ATIndex,'K2');
% 
% %% SD magnets
% ifam = ('SD');
% AO.(ifam).AT.ATType  = 'SEXT';
% AO.(ifam).AT.ATIndex = ATindx.SXD(:);
% AO.(ifam).Position   = s(AO.(ifam).AT.ATIndex);
% AO.(ifam).AT.ATParameterGroup = mkparamgroup(THERING,AO.(ifam).AT.ATIndex,'S');
% 
% % RF Cavity
% ifam = ('RF');
% AO.(ifam).AT.ATType = 'RF Cavity';
% AO.(ifam).AT.ATIndex = findcells(THERING,'Frequency')';
% AO.(ifam).Position   = s(AO.(ifam).AT.ATIndex);
setao(AO);
updateatindex;

if iscontrolroom
    switch2online;
else
    switch2sim;
end
