close all;
clear all;
file_name='20160204c_THOMX#009_fieldmap_200A';
%file_name='20160226b_THOMX#001_fieldmap_200A'
ext_in='.txt'
ext_out='_MP.txt'

%a = dlmread('20160204c_THOMX#009_fieldmap_200A.txt');
a = dlmread([file_name ext_in]);

%fid = fopen('lissage_lowess40.txt');
%nz=4001

a=a';
Svec=a(2:end,1)*1e3;
Xvec=a(1,2:end)*1e3;
Bz=abs(a(2:end,2:end));
%% 
%nn= find(Xvec==0)
ln=find(Svec==0)
%Bz0=Bz(nn,1:ln);
%Bz0=Bz(nn,:);

Svec=Svec(ln:end,1)
Bz=Bz(ln:end,:)
nn= find(Xvec==0)
Bz0=Bz(:,nn);
%%
%figure, plot(Svec(1:ln),Bz0);
figure, plot(Svec,Bz0);

%%

[S,X]=meshgrid(Svec(1):abs(Svec(1)-Svec(2)):Svec(end),Xvec(1):abs(Xvec(1)-Xvec(2)):Xvec(end));


figure    
     imagesc(Svec,Xvec,Bz');
      %imagesc(Bz(1,1:ln));
     shading interp;
     colorbar;
 	xlabel ('S(mm)');ylabel('X (mm)');
    
%%

VB=Bz'
rho0=352;
theta0=(21.5)*pi/180;

drho=-24:2:24;
%drho=0;
thetainf=pi/2;
thetasup=pi/2+theta0;
vartheta1=thetainf:0.001*pi/2:thetasup;
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
    
opera=dlmread('dipoleV33XC10.txt')
    

figure
    subplot(2,2,1); plot(DL_0,HO(:,4),'k.-')
    hold on, plot(opera(:,1),opera(:,2),'r')
    title('quadru(m^-1)')
    subplot(2,2,2); plot(DL_0,HO(:,3),'k.-')
    hold on, plot(opera(:,1),opera(:,3),'r')
    title('sectu(m^-2)')
    subplot(2,2,3); plot(DL_0,HO(:,2),'k.-')
    hold on, plot(opera(:,1),opera(:,4),'r')
    title('octu(m^-3)')
    subplot(2,2,4); plot(DL_0,HO(:,1),'k.-')
    hold on, plot(opera(:,1),opera(:,5),'r')
    title('deca(m^-4)')
V=[DL_0' HO(:,4) HO(:,3) HO(:,2) HO(:,1)]
%dlmwrite([file_name ext_out],V)

figure(4),
plot(DL_0,profilBz) 

%% 
[val, ni1]=min(abs(DL_0-69))
[val, ni2]=min(abs(DL_0-103.5))
[val, ni3]=min(abs(DL_0-138))
[val, ni4]=min(abs(DL_0-170))
[val, ni5]=min(abs(DL_0-210))
[val, ni6]=min(abs(DL_0-250))

ni=[1 126 189 251 315 385 478];%0.001
ni=[1 14 20 26 33 40 49];%0.01
ni=[ni1 ni2 ni3 ni4 ni5 ni6]
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
  sum(HL(:,4))
  
  
  
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
    %Refpart1=griddata(S,X,VB,-rho0*cos(vartheta1),rho0*sin(vartheta1)-rho0);
    Refpart1=interp2(S,X,VB,-rho0*cos(vartheta1),rho0*sin(vartheta1)-rho0,'spline');
    Refpart_0=[Refpart1 Refpart2];
    DL0=-(Ssup+sqrt((Xini-X0).*(Xini-X0)+(Sinf-S0).*(Sinf-S0)));
    DL2=(Sinf+sqrt((Xini-X2).*(Xini-X2)+(Sinf-S2).*(Sinf-S2)));
    DL1=rho0*(vartheta1-pi/2);
    DL_0=[DL1 DL2];
    Bx0=trapz(DL_0,Refpart_0);
    Lmag=abs(Bx0/max(abs(Refpart_0)))*1e-3

%% 
  
  
  
%%   
dx0=0.020;
DB_quadru=sum(HL(:,4))*2/Lmag*(rho0*1e-3)*dx0
DB_sextu=sum(HL(:,1))*2/Lmag*(rho0*1e-3)*dx0^2
DB_octu=sum(HL(:,2))*2/Lmag*(rho0*1e-3)*dx0^3
DB_deca=sum(HL(:,3))*2/Lmag*(rho0*1e-3)*dx0^4

