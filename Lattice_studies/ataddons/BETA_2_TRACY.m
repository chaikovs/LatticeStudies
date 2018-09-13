function BETA_2_TRACY
% Get quad & sextupole list from beta file "name"


% dir='/home/sources/physmach/loulergue/work/matlab/SOLEIL/Upgrade/Lattice-save/';
% input_file =[dir 'SOLEIL_U_v4.str'];
dir='/home/sources/physmach/loulergue/work/matlab/SOLEIL/Upgrade-collab/Lattice-16fold/';
input_file =[dir 'centre_opt1_7BA_type1_fit9_atmatch3.str'];
    
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

% Drift
L=[];
for i=1:nsdfamily
    N{i}=list_sd{i,1};
    L=[L,str2num(list_sd{i,3})];
end
fprintf('\n')
for i=1:nsdfamily
    fprintf('%5s  : drift, L=%14.6d ; \n',N{i},L(i))  
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
    fprintf('%5sa  : quadrupole, L=%14.6d/2,K=%14.6d*dgsurgL, FF1=quadfringe, FF2=0,FFscaling = scale, method=intmeth,N=Nq; \n',N{i},L(i),K(i))  
end
fprintf('\n')
for i=1:nqpfamily
    fprintf('%5sb  : quadrupole, L=%14.6d/2,K=%14.6d*dgsurgL,  FF1=0, FF2=quadfringe,FFscaling = scale, method=intmeth,N=Nq; \n',N{i},L(i),K(i))  
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
    fprintf('%5s : sextupole, L=%14.6d, K =%14.6d*coef, N=NqSx, FF1=sextfringe, FF2=sextfringe, method=method4sextu; \n',N{i},L(i),K(i))
end

% Dipole
L=[];T=[];K=[];
for i=1:ndifamily
    N{i}=list_di{i,1};
    L=[L,str2num(list_di{i,3})*str2num(list_di{i,4})];
    T=[T,str2num(list_di{i,3})*180/pi];
    K=[K,-str2num(list_di{i,5})/(str2num(list_di{i,4})^2)];
end

fprintf('\n')
for i=1:ndifamily
    fprintf('%5s  : bending, L=%14.6d, T=%14.6d, T1=0E0, T2=0E0, K=%14.6d, method=intmeth, gap=tracy_gap, N=4; \n',N{i},L(i),T(i),K(i))  
end

