function lattice_beta2AT
% Get quad & sextupole list from beta file "name"
% print into AT input file : copy and past 
% different elem naming ...
%

clear
% input and output beta files
%dir='/home/loulergue/work/structure/anneau/nanoscopium/optic-betax11m_SDL01_09/';
%input_file =[dir 'solnano_chic_155_229_14m_5m_11m_11m_realSX_1.2_2.0_test2.str'];
%input_file =[dir 'solnano_chic_155_229_14m_5m_11m_11m_realSX_1.2_1.8_test2_3.7nm.str'];

% dir='/home/loulergue/work/structure/anneau/nanoscopium/optic_betax11m_SDL01_09_13/';
% input_file =[dir 'solnano_chic_155_229_14m_5m_11m_11m_realSX_1.2_2.0_test25µrad-fit2-chro1.str'];

% dir='/home/sources/physmach/loulergue/work/structure/anneau/nanoscopium/optic_betaz1m_betax11m_SDL01_09_13/';
% input_file =[dir 'OPA/solnano_chic_160_225_12m_5m_12m_12m_realSX_1.2_2.0_bz1.0m_opt2_opa6_test6_limQ11.str'];

% dir ='/home/sources/physmach/loulergue/work/matlab/Simulation/beta_structure/';
% input_file =[dir 'Upgrade/6BA-7BA-AB/Full-ring-2/6BA-DLS-scaled-2-sx2-fit1-optim1.str'];

dir ='/home/sources/physmach/loulergue/work/matlab/Simulation/beta_structure/';
input_file =[dir 'Upgrade/test-offmomentum/6BA-DLS-3Q10-fit.str'];
%input_file =[dir 'Upgrade/test-offmomentum/7BA-ESRF-3Q2-fit.str'];

fprintf('##########################\n')
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lecture input file
nfile=size(input_file);nfile=nfile(2);

clear list_qp  list_sx
fichier=input_file;
fid = fopen(fichier,'r');
if fid==-1
    display('Bug file')
    return
else
    fprintf('Input file %s \n',fichier)
end
nline=0;
% stock lines in txt
while (~feof(fid))
    C = textscan(fid,'%s',1,'delimiter','\n');
    nline=nline+1;
    txt(nline)=C;
end
fclose(fid);
nline=nline-1;

% Get blocks
block_type={'*** VERSION ***', '*** AUTHOR ***','*** TITRE ***',...
    '*** LIST OF ELEMENTS ***','*** STRUCTURE ***','*** PERIOD ***',...
    '*** OPTION ***','*** BEAM-MATRIX ***','*** DISPERSION ***',...
    '*** PARTICLE TYPE ***','*** ENERGIE CINETIQUE (MeV) ***',...
    '*** EMITTANCE ***','*** PARAMETERS OF FIT ***','*** VARIABLES ***',...
    '*** CONSTRAINTS ***','*** ENDFILE ***'};
nblock=length(block_type);
for i=1:nline
    for j=1:nblock
        test=strcmp(txt{i},block_type(j));
        if test
            block_line(j)=i;
            %fprintf('block %32s      at line %d  \n',block_type{j},block_line(j))
        end
    end
end


% Get quad & sext family
% in block *** LIST OF ELEMENTS ***, n=4
nl1=block_line(4)+2;
nl2=block_line(5)-1;
nsdfamily=0;
nqpfamily=0;
nsxfamily=0;
ndifamily=0;
for i=nl1:nl2
    ligne=txt{i};
    list=strread(char(ligne),'%s');
    if strcmp(list{2},'SD')  % Drift
        nsdfamily=nsdfamily+1;
        list_sd(nsdfamily,:)=list;
    end
    if (length(list)>1) && (strcmp(list{2},'QP')) ;
        nqpfamily=nqpfamily+1;
        list_qp(nqpfamily,:)=list;
        txt{i}=[list{1} ,'        ',list{2},'  ',list{3},'  ',list{4},'  ',list{5}];
    end
    if (length(list)>1) && (strcmp(list{2},'SX')) ;
        nsxfamily=nsxfamily+1;
        list_sx(nsxfamily,:)=list;
    end
    if strcmp(list{2},'DI')  % dipole
        ndifamily=ndifamily+1;
        list_di(ndifamily,:)=list;
    end

end

%Drift
L=[];
for i=1:nsdfamily
    N{i}=list_sd{i,1};
    L=[L,str2num(list_sd{i,3})];
end

fprintf('\n')
for i=1:nsdfamily
    fprintf('%5s = drift(''%3s'' , %14.6d, ''DriftPass''); \n',N{i},N{i},L(i))  
end

% Quadrupoles
L=[];K=[];
for i=1:nqpfamily
    N{i}=list_qp{i,1};
    L=[L,str2num(list_qp{i,3})];
    K=[K,str2num(list_qp{i,4})];
end

fprintf('\n')
for i=1:nqpfamily
    fprintf('%5s = quadrupole(''%5s'' , %14.6d, %14.6d , QPassMethod); \n',N{i},N{i},L(i),K(i))  
end


% Sextupoles
L=[];K=[];
for i=1:nsxfamily
    N{i}=list_sx{i,1};
    L=[L,str2num(list_sx{i,3})];
    K=[K,str2num(list_sx{i,4})];
end

fprintf('\n')
for i=1:nsxfamily
    fprintf('%5s =  sextupole( ''%3s'' , %14.6d,  %14.6d, SPassMethod); \n',N{i},N{i},L(i),K(i))
end

% dipole
L=[];T=[];K=[];
for i=1:ndifamily
    N{i}=list_di{i,1};
    L=[L,str2num(list_di{i,3})*str2num(list_di{i,4})];
    T=[T,str2num(list_di{i,3})];
    K=[K,-str2num(list_di{i,5})/(str2num(list_di{i,4})^2)];
end
%BEND  =  rbend2('BEND', 1.05243, theta, thetae, thetas, K,fullgap,'BndMPoleSymplectic4Pass');
fprintf('\n')
for i=1:ndifamily
    fprintf('%5s = rbend2(''%3s'' , %14.6d, %14.6d, %3.1d, %3.1d, %14.6d, %3.1d, ''BndMPoleSymplectic4Pass''); \n',N{i},N{i},L(i),T(i),0,0,K(i),0); 
end
