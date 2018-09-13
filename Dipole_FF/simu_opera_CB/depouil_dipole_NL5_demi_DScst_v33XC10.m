close all;
clear all;
fid = fopen('dipole33_XC10.table');
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
theta0=(22.5)*pi/180;

%%

figure    
     %imagesc(Svec,Xvec,VB);
     imagesc(Svec(250:750),Xvec,VB(:,250:750));
     shading interp;
     colorbar;
 	xlabel ('S(mm)');ylabel('X (mm)');
%%
drho=-24:2:24;
%drho=0;
thetainf=pi/2;
thetasup=pi/2+theta0;
vartheta1=thetainf:0.01*pi/2:thetasup;
rhoL=rho0+drho;
Xini=rhoL*sin(pi/2-theta0)-rho0;
Sinf=-rhoL*cos(pi/2-theta0);
Ssup=rhoL*cos(pi/2-theta0);

for j=1:length(vartheta1)
    X1=rhoL*sin(vartheta1(j))-rho0;
    S1=-rhoL*cos(vartheta1(j));
    Refpart=interp2(S,X,VB,S1,X1,'spline');
    if(j==1)
        profilBz1=Refpart';
    else
        profilBz1=[profilBz1 Refpart'];
    end
end
DL=abs(-rho0*(cos(vartheta1(1))-cos(vartheta1(2))));
nk=int16((260-rho0*cos(pi/2-theta0))/DL);
for k=1:nk
    S2=Ssup+double(k)*DL;
    X2=-S2.*(tan(theta0))+Xini-Sinf*tan(theta0);
    Refpart2=interp2(S,X,VB,S2,X2,'spline');
    if(k==1)
        profilBz2=Refpart2';
    else
        profilBz2=[profilBz2 Refpart2'];
    end
end

profilBz1;
profilBz2;
profilBz=[profilBz1 profilBz2];
B0=interp2(S,X,VB,0,0,'spline');
for i=1:length(vartheta1)+nk
    DB(:,i)=(profilBz(:,i)-profilBz(length(drho)/2+0.5,i))./B0;
    HO(i,1:5)=polyfit(drho*1e-3,DB(:,i)',4);
end
Xini=rho0*sin(pi/2-theta0)-rho0;
Sinf=-rho0*cos(pi/2-theta0);
Ssup=rho0*cos(pi/2-theta0);
S2=(Ssup)+DL:DL:(Ssup)+double(nk)*DL;
X2=-S2.*(tan(theta0))+Xini-Sinf*tan(theta0);
DL2=(Sinf+sqrt((Xini-X2).*(Xini-X2)+(Sinf-S2).*(Sinf-S2)));
DL1=rho0*(vartheta1-pi/2);
DL_0=[DL1 DL2];
    
    

figure
    subplot(2,2,1); plot(DL_0,HO(:,4),'k.-')
    title('quadru(m^-1)')
    subplot(2,2,2); plot(DL_0,HO(:,3),'k.-')
    title('sectu(m^-2)')
    subplot(2,2,3); plot(DL_0,HO(:,2),'k.-')
    title('octu(m^-3)')
    subplot(2,2,4); plot(DL_0,HO(:,1),'k.-')
    title('deca(m^-4)')
V=[DL_0' HO(:,4) HO(:,3) HO(:,2) HO(:,1)]
dlmwrite('dipoleV33XC10.txt',V)


%% figure IPAC Cynthia 2016

figure
    subplot(2,2,1); plot(DL_0(1:end-1)',HO(1:end-1,4),'k.-')
    ylabel('K_2/2(m^{-1})'), xlabel ('s(mm)')
    subplot(2,2,2); plot(DL_0(1:end-1)',HO(1:end-1,3),'k.-')
    ylabel('K_3/6(m^{-2})'), xlabel ('s(mm)')
    subplot(2,2,3); plot(DL_0(1:end-1)',HO(1:end-1,2),'k.-')
    ylabel('K_4/24(m^{-3})'), xlabel ('s(mm)')
    subplot(2,2,4); plot(DL_0(1:end-1)',HO(1:end-1,1),'k.-')
    ylabel('K_5/120(m^{-4})'), xlabel ('s(mm)')
    
set(gcf,'Color','white')
%set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)

%% figure IEEE Cynthia 2017

figure
    subplot(2,2,1); plot([-DL_0(end-1:-1:2) DL_0(1:end-1)]',[HO(end-1:-1:2,4)' HO(1:end-1,4)'],'k.-')
    ylabel('K_2/2(m^{-1})'), xlabel ('s(mm)')
    subplot(2,2,2); plot([-DL_0(end-1:-1:2) DL_0(1:end-1)]',[HO(end-1:-1:2,3)' HO(1:end-1,3)'],'k.-')
    ylabel('K_3/6(m^{-2})'), xlabel ('s(mm)')
    subplot(2,2,3); plot([-DL_0(end-1:-1:2) DL_0(1:end-1)]',[HO(end-1:-1:2,2)' HO(1:end-1,2)'],'k.-')
    ylabel('K_4/24(m^{-3})'), xlabel ('s(mm)')
    subplot(2,2,4); plot([-DL_0(end-1:-1:2) DL_0(1:end-1)]',[HO(end-1:-1:2,1)' HO(1:end-1,1)'],'k.-')
    ylabel('K_5/120(m^{-4})'), xlabel ('s(mm)')

Vtot=[[-DL_0(end-1:-1:2) DL_0(1:end-1)]' [HO(end-1:-1:2,4)' HO(1:end-1,4)']' [HO(end-1:-1:2,3)' HO(1:end-1,3)']' [HO(end-1:-1:2,2)' HO(1:end-1,2)']' [HO(end-1:-1:2,1)' HO(1:end-1,1)']']
dlmwrite('dipoleV33XC10_tot.txt',Vtot)
   
set(gcf,'Color','white')
%set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)

%%
ni=[1 126 189 251 315 385 478];%0.001
ni=[1 14 20 26 33 40 49];%0.01
    for i=1:length(ni)-1
        HL(i,1)=trapz(DL_0(1,ni(i):ni(i+1)), HO(ni(i):ni(i+1),3))/rho0;
        HL(i,2)=trapz(DL_0(1,ni(i):ni(i+1)), HO(ni(i):ni(i+1),2))/rho0;
        HL(i,3)=trapz(DL_0(1,ni(i):ni(i+1)), HO(ni(i):ni(i+1),1))/rho0;
        HL(i,4)=trapz(DL_0(1,ni(i):ni(i+1)), HO(ni(i):ni(i+1),4))/rho0;

    end
    format long
    HL
  sum(HL(:,1))
  sum(HL(:,2))
  sum(HL(:,3))
  
        HLtot(1)=trapz(DL_0, HO(:,3))/rho0;
        HLtot(2)=trapz(DL_0, HO(:,2))/rho0;
        HLtot(3)=trapz(DL_0, HO(:,1))/rho0;
        HLtot(4)=trapz(DL_0, HO(:,4))/rho0;
%% Calcul de Lmag
    Xini=rho0*sin(pi/2-theta0)-rho0;
    Sinf=-rho0*cos(pi/2-theta0);
    Ssup=rho0*cos(pi/2-theta0);
    S0=-250:1:(Sinf-0.5);
    S2=(Ssup+0.5):1:250;
    X0=S0.*(tan(theta0))+Xini-Sinf*tan(theta0);
    X2=-S2.*(tan(theta0))+Xini-Sinf*tan(theta0);
    Stot=[-rho0*cos(vartheta1) S2];
    Xtot=[rho0*sin(vartheta1)-rho0 X2];
    Refpart2=interp2(S,X,VB,S2,X2,'spline');
    Refpart1=interp2(S,X,VB,-rho0*cos(vartheta1),rho0*sin(vartheta1)-rho0,'spline');
    Refpart_0=[Refpart1 Refpart2];
    DL0=-(Ssup+sqrt((Xini-X0).*(Xini-X0)+(Sinf-S0).*(Sinf-S0)));
    DL2=(Sinf+sqrt((Xini-X2).*(Xini-X2)+(Sinf-S2).*(Sinf-S2)));
    DL1=rho0*(vartheta1-pi/2);
    DL_0=[DL1 DL2];
    Bx0=trapz(DL_0,Refpart_0);
    Lmag=abs(Bx0/max(abs(Refpart_0)))*1e-3

  
  
%%   
dx0=0.020;
DB_quadru=sum(HL(:,4))*2/Lmag*(rho0*1e-3)*dx0
DB_sextu=sum(HL(:,1))*2/Lmag*(rho0*1e-3)*dx0^2
DB_octu=sum(HL(:,2))*2/Lmag*(rho0*1e-3)*dx0^3
DB_deca=sum(HL(:,3))*2/Lmag*(rho0*1e-3)*dx0^4



DB_quadru_tot=HLtot(4)*2/Lmag*(rho0*1e-3)*dx0
DB_sextu_tot=HLtot(1)*2/Lmag*(rho0*1e-3)*dx0^2
DB_octu_tot=HLtot(2)*2/Lmag*(rho0*1e-3)*dx0^3
DB_deca_tot=HLtot(3)*2/Lmag*(rho0*1e-3)*dx0^4


