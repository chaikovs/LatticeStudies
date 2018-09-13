close all;
clear all;
%fid = fopen('dipole33_XC10.table');
fid = fopen('dipole33_XC10_2016_200A.table')
for i=1:7
ligne=fgetl(fid);
end
a = fscanf(fid,'%f\t %f\t %f\t %f\n',[4 Inf]);
fclose(fid);
a=a';
Xvec=a(:,1);
Svec=a(:,3);
Bz=a(:,4);
D=0.8;
VB=Bz(1:1001,1)';
for i=2:201
B_map=Bz((i-1)*1001:(i-1)*1001+1000,1);
i
length(B_map);
VB=[VB ;B_map'];
end
[S,X]=meshgrid(-500:1:500,-100:1:100);

rho0=352;
theta0=(21.5)*pi/180;

figure    
     imagesc(Svec,Xvec,VB);
     shading interp;
     colorbar;
 	xlabel ('S(mm)');ylabel('X (mm)');

rho0=352;
%theta0=(22.5)*pi/180;
theta0=(21.5)*pi/180;

drho=-24:2:24;
drho=0;
thetainf=pi/2;
thetasup=pi/2+theta0;
vartheta1=thetainf:0.01*pi/2:thetasup;
rhoL=rho0+drho;

%% curved trajectory
for j=1:length(vartheta1)
    X1(j)=rhoL*sin(vartheta1(j))-rho0;
    S1(j)=-rhoL*cos(vartheta1(j));
    Refpart_cpos(j)=interp2(S,X,VB,S1(j),X1(j)); % positive s part
    Refpart_cneg(j)=interp2(-S,X,VB,S1(j),X1(j)); %negative s part

%     if(j==1)
%         profilBz_cpos=Refpart_cpos';
%         profilBz_cneg=Refpart_cneg';
%        
%     else
%         profilBz_cpos=[profilBz_cpos Refpart_cpos'];
%         profilBz_cneg=[profilBz_cneg Refpart_cneg'];
%     end
end


DL=abs(-rho0*(cos(vartheta1(1))-cos(vartheta1(2))));
nk=int16((Svec(end)-rho0*cos(pi/2-theta0))/DL);
Slin=-rho0*cos(vartheta1(end)):DL:300
for k=1:length(Slin)
    S2(k)=Slin(k);
    X2(k)=X1(end)-(Slin(k)-S1(end)).*(tan(theta0));
    Refpart_lpos(k)=interp2(S,X,VB,S2(k),X2(k),'spline');
    Refpart_lneg(k)=interp2(-S,X,VB,S2(k),X2(k),'spline');
%     if(k==1)
%         profilBz_lpos=Refpart_lpos';
%         profilBz_lneg=Refpart_lneg';
% 
%     else
%         profilBz_lpos=[profilBz_lpos Refpart_lpos'];
%         profilBz_lneg=[profilBz_lneg Refpart_lneg'];
%     end
end

%profilBz=[profilBz_lneg(:,end:-1:1) profilBz_cneg(:,end:-1:1) profilBz_cpos profilBz_lpos];


DL_0=[-S2(end:-1:1) -S1(end:-1:1) S1 S2];
%figure, plot(DL_0,profilBz)    
figure(4), plot([-S2(end:-1:1) -S1(end:-1:1) S1 S2],[Refpart_lneg(:,end:-1:1) Refpart_cneg(:,end:-1:1) Refpart_cpos Refpart_lpos] )
dlmwrite('profilBz_operaV33XC10_215deg.txt',[[-S2(end:-1:1) -S1(end:-1:1) S1 S2]' [Refpart_lneg(:,end:-1:1) Refpart_cneg(:,end:-1:1) Refpart_cpos Refpart_lpos]'],'delimiter','\t')
%% lecture fichier opera

