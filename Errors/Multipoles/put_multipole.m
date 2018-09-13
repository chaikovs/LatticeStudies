function put_multipole
% Get beta file "name"
%     must be with full length quad and sext
%     No already multipole in structure, otherwise added
% Add multipole list
% Split quad and add multipole in the middle
% Add sext multipoleust after (no split yet)
% Create new file "name_multipole"

clear
% directories
dir='/Users/ichaikov/Documents/MATLAB/thomx-mml/machine/THOMX/StorageRing/Errors/Multipoles/';
% input and output beta files
%input_file  =[dir 'TDR_0.17_0.64_r56_0.2_sx_Dff_Dsx_Doc_D10_v31.str'];
%output_file =[dir 'TDR_0.17_0.64_r56_0.2_sx_Dff_Dsx_Doc_D10_v31_mult_test.str'];
input_file  =[dir 'TDR_0.17_0.64_r56_0.2_sx_Dff41.2_FF_chro00_SXfit1_Qff_DipM.str'];
output_file =[dir 'TDR_0.17_0.64_r56_0.2_sx_Dff41.2_FF_chro00_SXfit1_Qff_DipQuadSextM.str'];



% input multipole files
qp_file='qp_multipoles_mean_meas.txt';
sx_file='sx_multipoles.txt';


fprintf('##########################\n')
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lecture input file
fid = fopen(input_file,'r');
if fid==-1
    display('Bug file')
    return
else
    fprintf('Input file %s \n',input_file)
end
nline=0;
% stock lines in txt
while (~feof(fid)) 
    C = textscan(fid,'%s',1,'delimiter','\n');
    nline=nline+1;
    txt(nline)=C;
end
fclose(fid);
nline=nline%;-1;

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
          fprintf('block %32s      at line %d  \n',block_type{j},block_line(j))
       end
    end
end
fprintf('\n')

% Get quad & sext family & divide length by 2
% in block *** LIST OF ELEMENTS ***, n=4
nl1=block_line(4)+2;
nl2=block_line(5)-1;
nqpfamily=0;
nsxfamily=0;
for i=nl1:nl2
    ligne=txt{i};
    list=strread(char(ligne),'%s');
    if (length(list)>1) && (strcmp(list{2},'QP')) ;
        nqpfamily=nqpfamily+1;
        list_qp(nqpfamily,:)=list;
        L=str2num(list{3});L=L/2;
        list{3}=num2str(L);
        txt{i}=[list{1} ,'        ',list{2},'  ',list{3},'  ',list{4},'  ',list{5}];
    end
     if (length(list)>1) && (strcmp(list{2},'SX')) ;
        nsxfamily=nsxfamily+1;
        list_sx(nsxfamily,:)=list;
     end
end
fprintf('Found %d quad families    \n',nqpfamily)
fprintf('Found %d sext families    \n',nsxfamily)
fprintf('\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate multipoles from list for each families
fid = fopen(qp_file,'r');
nqp_multipoles=0;
while (~feof(fid)) 
    C = textscan(fid,'%s',1,'delimiter','\n');
    nqp_multipoles=nqp_multipoles+1;
    qp_multipoles_list(nqp_multipoles)=C;
    if length( char(C{1}))==0
        nqp_multipoles=nqp_multipoles-1;
        break
    end
end
fclose(fid);

fprintf('Found %d quad multipoles   \n',nqp_multipoles)
fprintf('   Skew   order   Strength \n')
for i=1:nqp_multipoles
    ligne=qp_multipoles_list{i};
    list=strread(char(ligne),'%s');
    type       =char(list(1));
    skew(i)    =type(1);
    ordre(i)   =str2num(char(type(2:length(type))));
    strength(i)=str2num(char(list(2)));
    fprintf('    %s     %3i     %6.2f    \n',skew(i),ordre(i),strength(i)*1e4)
end
fprintf('\n')

% Make QP multipole for beta
x0=0.02;  % mesure ??? 20 mm
k=1;
for i=1:nqpfamily
    qp_name    =list_qp(i,1);
    qp_length  =str2num(char(list_qp(i,3)));
    qp_strength=str2num(char(list_qp(i,4)));
    for j=1:nqp_multipoles
        type='LD';
        if skew(j)=='a' ; type='LT'; end;
        name=strcat(type, num2str(ordre(j)) ,qp_name);
        npole=2*ordre(j);
        L_strength=strength(j)*qp_length*qp_strength/power(x0,ordre(j)-2);% calcul force quad
        var  =0.00;
        txt_qp(i,j,1:5)=[name type num2str(L_strength) num2str(npole)  num2str(var)];
        %fprintf('%s  %s  %s  %s  \n',char(txt_qp(i,j,1)),char(txt_qp(i,j,2)),char(txt_qp(i,j,3)),char(txt_qp(i,j,4)))
        k=k+1;
    end
end


fid = fopen(sx_file,'r');
nsx_multipoles=0;
while (~feof(fid)) 
    C = textscan(fid,'%s',1,'delimiter','\n');
    nsx_multipoles=nsx_multipoles+1;
    sx_multipoles_list(nsx_multipoles)=C;
    if length( char(C{1}))==0
        nsx_multipoles=nsx_multipoles-1;
    break
    end

end
fclose(fid);

fprintf('Found %d sext multipoles   \n',nsx_multipoles)
fprintf('   Skew   order   Strength   \n')
for i=1:nsx_multipoles
    ligne=sx_multipoles_list{i};
    list=strread(char(ligne),'%s');
    type       =char(list(1));
    skew(i)    =type(1);
    ordre(i)   =str2num(char(type(2:length(type))));
    strength(i)=str2num(char(list(2)));
    fprintf('    %s     %3i     %6.2f  \n',skew(i),ordre(i),strength(i)*1e4)
end
fprintf('\n')


% Make SX multipole for beta
x0=0.020;  % mesure ??? 20 mm
k=1;
for i=1:nsxfamily
    sx_name    =list_sx(i,1);
    sx_length  =str2num(char(list_sx(i,3)));
    sx_strength=str2num(char(list_sx(i,4)));
    for j=1:nsx_multipoles
        type='LD';
        if skew(j)=='a' ; type='LT'; end;
        name=strcat(type, num2str(ordre(j)) ,sx_name);
        npole=2*ordre(j);
        L_strength=strength(j)*sx_length*sx_strength/power(x0,ordre(j)-3); % calcul force sextu
        var  =0.00;
        txt_sx(i,j,1:5)=[name type num2str(L_strength) num2str(npole)  num2str(var)];
        %fprintf('%s  %s  %s  %s  \n',char(txt_sx(i,j,1)),char(txt_sx(i,j,2)),char(txt_sx(i,j,3)),char(txt_sx(i,j,4)))
        k=k+1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write new file
block_line(nblock+1)=block_line(nblock)+1;
fid = fopen(output_file,'w');
for i=1:nblock
    j1=block_line(i);
    j2=block_line(i+1)-1;
    if i~=5
        for k=j1:j2
            fprintf(fid,' %s \n',char(txt{k}));
        end
        % on rajoute les multipoles sur list of element
        if i==4
            for i=1:nqpfamily
                for j=1:nqp_multipoles
                    fprintf(fid,' %s    %s    %s    %s    %s\n',...
                        char(txt_qp(i,j,1)),char(txt_qp(i,j,2)),...
                        char(txt_qp(i,j,3)),char(txt_qp(i,j,4)),char(txt_qp(i,j,5)));
                end
            end
            for i=1:nsxfamily
                for j=1:nsx_multipoles
                    fprintf(fid,' %s    %s    %s    %s    %s\n',...
                        char(txt_sx(i,j,1)),char(txt_sx(i,j,2)),...
                        char(txt_sx(i,j,3)),char(txt_sx(i,j,4)),char(txt_sx(i,j,5)));
                end
            end
        end
    else
        % on rajoute les multipoles sur la structure
        nnqp=0; nnsx=0;
        fprintf(fid,' %s \n',char(txt{j1}));
        fprintf(fid,' %s \n',char(txt{j1+1}));
        for k=j1+2:j2
            ligne=txt{k};
            list0=strread(char(ligne),'%s');list=[];
            ll=length(list0);
            
            for l=1:ll
                elem=list0{l};
                flag=0;
                for nf=1:nqpfamily % recherche QP
                    if strcmp(elem,list_qp{nf,1})
                        nnqp=nnqp+1;
                        nf1=nf;
                        flag=1;
                    end
                end
                for nf=1:nsxfamily % recherche SX
                    if strcmp(elem,list_sx{nf,1})
                        nnsx=nnsx+1;
                        nf2=nf;
                        flag=2;
                    end
                end
                list=[list ,'  ', elem];
                if flag==1
                    fprintf(fid,' %s \n',list);list=[];
                    %lentille
                    for jj=1:nqp_multipoles
                        lentille=char(txt_qp(nf1,jj,1));
                        list=[list ,'  ', lentille ];
                    end
                    list=[list ,'  ', elem];
                end
                if flag==2
                    fprintf(fid,' %s \n',list);list=[];
                    %lentille
                    for jj=1:nsx_multipoles
                        lentille=char(txt_sx(nf2,jj,1));
                        list=[list ,'  ', lentille ];
                    end
                end
                
                
            end
            % coupe la list avant ecriture
            %length(list);
            
            fprintf(fid,' %s \n',list);
            
        end
    end
end
fprintf('Found %d quad    \n',nnqp)
fprintf('Found %d sext    \n',nnsx)
fprintf('\n')

fclose(fid);
fprintf('output file %s  amended with multipoles\n',output_file)



% fprintf('Found %d quad multipoles   \n',nqp_multipoles)
% fprintf('   Skew   order   Strength_short   Strength_long  \n')
% for i=1:nqp_multipoles
%     ligne=qp_multipoles_list{i};
%     list=strread(char(ligne),'%s');
%     type       =char(list(1));
%     skew(i)    =type(1);
%     ordre(i)   =str2num(char(type(2:length(type))));
%     strength_short(i)=str2num(char(list(2)));
%     strength_long(i)=str2num(char(list(3)));
%     fprintf('    %s     %3i     %6.2f          %6.2f \n',skew(i),ordre(i),strength_short(i)*1e4,strength_long(i)*1e4)
% end
% fprintf('\n')
% 
% % Make QP multipole for beta
% L_short=0.32;
% L_long =0.46;
% x0=0.030;  % mesure ??? 30 mm
% k=1;
% for i=1:nqpfamily
%     qp_name    =list_qp(i,1);
%     qp_length  =str2num(char(list_qp(i,3)));
%     qp_strength=str2num(char(list_qp(i,4)));
%     for j=1:nqp_multipoles
%         type='LD';
%         if skew(j)=='a' ; type='LT'; end;
%         name=strcat(type, num2str(ordre(j)) ,qp_name);
%         npole=2*ordre(j);
%         strength=strength_short(j);
%         if qp_length==L_long ; strength=strength_long(j); end ;
%         L_strength=strength*qp_length*qp_strength/power(x0,ordre(j)-2);% calcul force quad
%         var  =0.00;
%         txt_qp(i,j,1:5)=[name type num2str(L_strength) num2str(npole)  num2str(var)];
%         %fprintf('%s  %s  %s  %s  \n',char(txt_qp(i,j,1)),char(txt_qp(i,j,2)),char(txt_qp(i,j,3)),char(txt_qp(i,j,4)))
%         k=k+1;
%     end
% end
% 
% 
