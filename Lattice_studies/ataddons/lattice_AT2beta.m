% basic translator from AT to beta file
%
%clear
% cd base-ryutaro
% %RING=lat_7BA_Ryutaro_2;
% RING=lat_6BA_Ryutaro_2;
% cd ..
%
len=length(RING);
fprintf('lattice length : %4i \n', len)
fprintf('\n')

%
maille={};
for i=1:len
    class=RING{i}.Class;
    if strcmp(class,'Drift')    
        fprintf('%s    %s  %d \n',RING{i}.FamName,'SD',RING{i}.Length)
    end
    if strcmp(class,'Quadrupole');
        fprintf('%s    %s  %d  %d\n',RING{i}.FamName,'QP',RING{i}.Length,RING{i}.K)
    end
    if strcmp(class,'Bend') 
        R=RING{i}.Length/RING{i}.BendingAngle;
        ind=abs(RING{i}.K*R*R);
        fprintf('%s    %s  %d  %d  %d\n',RING{i}.FamName,'DI',RING{i}.BendingAngle,R,ind)
    end
    
    maille{i}=RING{i}.FamName;
    
end
fprintf('\n')

for k=1:6:len
    fprintf('%s  %s  %s  %s  %s  %s \n',maille{k:k+5})
end


