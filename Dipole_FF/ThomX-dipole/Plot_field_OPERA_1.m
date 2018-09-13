%function  Plot_field_OPERA
%PLOT_FIELD Summary of this function goes here
%   Detailed explanation goes here
%   Maille ThomX
% Hard hedge traj

%Dipole 
radius  =0.352; %m
teta    =pi/8; % half field deviation
tetapole=pi/9; % half yoke deviation
gap     =0.042 ;

fprintf(' \n')
fprintf('######################## \n')
fprintf('Some plot from OPERA mTm \n')

%[Xb ,Sb , Bz] = getOPERAfield2D('field/champ_dipole_OPERA.table');
%[Xb ,Sb , Bz] = getOPERAfield2D('field/champ_dipole_OPERA_v6.table');
%[Xb ,Sb , Bz] = getOPERAfield2D('field/champ_dipole_OPERA_v33.table');
%[Xb ,Sb , Bz] = getOPERAfield2D('field/champ_dipole_OPERA_v33bis.table');

%[Xb ,Sb , Bz] = getOPERAfield2D('field/dipole33_XC10.table');
[Xb,Sb,Bz] = getMEASUREDfield2D('field/20160208a_THOMX#009_fieldmap_160A.xlsx');
Bz=Bz*0.996; % *0.996 pour ~50 MeV 
%[Xb ,Sb , Bz] = getOPERAfield2D('field/dipole-perfect_XC10.table');
%[Xb ,Sb , Bz] = getOPERAfield2D('field/dipole-defaut_XC10.table');
%[Xb,Sb,Bz] = getMEASUREDfield2D('field-meas/20160129a_THOMX#009_fieldmap_300A.xlsx');
bz0=interp2(Xb,Sb,Bz,0,0);
fprintf('   Peak field  %g  T \n',bz0 )

%
[ traj0 , ds, np  ] = make_traj(radius,0,teta);
bz   =interp2(Xb,Sb,Bz,traj0(1,:),traj0(2,:));
ds0=ds(1);
Lmag  =sum(bz.*ds)/bz(1); % equivalent half length

bz   =interp2(Xb,Sb,Bz,traj0(1,:)  ,traj0(2,:));
% bzp20=interp2(Xb,Sb,Bz,trajp20(1,:),trajp20(2,:));
% bzm20=interp2(Xb,Sb,Bz,trajm20(1,:),trajm20(2,:));
intbz0=sum(bz.*ds);
intbz00=trapz(cumsum(ds),bz)*1000;
fprintf('  Integrated field  %g  mTm \n',intbz00 *2)

% Get profil along path bzn
xx=(-24:2:24);
intbz=[];
bzn=[];
for dx=xx
    [ traj , ds ] = make_traj((radius+dx*1e-3),dx*1e-3,teta);
    bz=interp2(Xb,Sb,Bz,traj(1,:),traj(2,:));
    intbz=[intbz sum(bz.*ds)];
    bzn=[bzn ; bz.*ds/ds0];
    %bzn=[bzn ; bz];
end


% Get poly fit on integrated form
intbzr =(intbz-intbz0)/intbz0;
intbzrc=(intbz-intbz0 - bz0*tetapole*xx/1000)/intbz0;
p = polyfit(xx*1e-3,intbzrc,4);pint=p;


%%

f = polyval(pint,xx*1e-3);
figure(11)
set(gca,'Fontsize',18)
plot(xx,intbzr*100,'ob');hold on
plot(xx,intbzrc*100,'or');
plot(xx,f*100,'--r'); hold off
xlabel('dx (mm)')
ylabel('d intB / intB (%)')
legend('Direct int field profile','Corr int field profile', 'Fit')
grid on

% Get polyfit along path 
pn=[];
% for i=1:1:length(ds)
%     db=(bzn(:,i)-bzn(11,i))/bzn(11,1);
%     p = polyfit(xx*1e-3,db',4);
%     pn=[pn ; p];
% end

for i=1:1:length(ds)
    corr=0;% to remove natural dip H foc in matrix
    if i<(np+2) 
        corr=bzn(13,i)*xx'*1e-3*teta/(np+1)/ds0;
    end % before 11
    db=(bzn(:,i)-bzn(13,i)-corr)/bzn(13,1);
    p = polyfit(xx*1e-3,db',4);
    pn=[pn ; p];
end



% get polyfit sum
% sum bit diff from above : ds variation not included !
ps=sum(pn)*ds0/Lmag;

ff = polyval(ps,xx*1e-3);

%%
rx=2.4e-2; 

DBB_quad = (ps(:,4))*(rx);
DBB_sext = (ps(:,3))*(rx)^2;
DBB_oct = (ps(:,2))*(rx)^3;
DBB_deca = (ps(:,1))*(rx)^4;

Bint_quad = DBB_quad*bzn(13,1)*2*Lmag;
Bint_sext = DBB_sext*bzn(13,1)*2*Lmag;
Bint_oct = DBB_oct*bzn(13,1)*2*Lmag;
Bint_deca = DBB_deca*bzn(13,1)*2*Lmag;

fprintf('   Int B2 quad [T m]  %10.5E \n',Bint_quad)
fprintf('   Int B3 sext [T m]%10.5E  \n',Bint_sext)
fprintf('   Int B4 octo [T m] %10.5E  \n',Bint_oct)
fprintf('   Int B5 deca [T m] %10.5E  \n',Bint_deca)

B_quad = pn(:,4)*(rx)*bzn(13,1);
B_sext = pn(:,3)*(rx)^2*bzn(13,1);
B_oct = pn(:,2)*(rx)^3*bzn(13,1);
B_deca = pn(:,1)*(rx)^4*bzn(13,1);
%%

[x,Bint,dBB] = getMEASUREDMultipoles('field/20160208a_THOMX#009_fieldmap_160A.xlsx');

figure(1)
set(gca,'Fontsize',18)
plot(xx,polyval(ps,xx*1e-3),'ob');hold on 
plot(xx,polyval(pint,xx*1e-3),'or');
plot(x,dBB,'ok');
legend('polyfit sum','polyfit on integrated form', 'ALBA results')
xlabel('dx (mm)')
ylabel('dB/B ')
hold off 
grid on

%%

figure
set(gca,'Fontsize',18)
plot(xx,Bint,'*-');
xlabel('dx (mm)')
ylabel('Bint [T m] ')

%%
% [s_alba,b2_alba,b3_alba,b4_alba,b5_alba] = getMEASUREDMultipolesSep('field/20160208a_THOMX#009_fieldmap_160A.xlsx');
% 
% % plot grad
% figure(22)
% plot(traj0(2,:),B_quad,'-r');hold on
% plot(s_alba,b2_alba,'-b')
% hold off
% xlabel('S (mm)')
% ylabel('Gradient B_2 [T@Ref]')
% grid on
% 
% % plot sext
% figure(33)
% plot(traj0(2,:),B_sext,'-r');hold on
% plot(s_alba,b3_alba,'-b')
% xlabel('S (mm)')
% ylabel('Sext  B_3 [T@Ref]')
% grid on
% 
% 
% % plot oct
% figure(44)
% plot(traj0(2,:),B_oct,'-r');hold on
% plot(s_alba,b4_alba,'-b')
% xlabel('S (mm)')
% ylabel('Oct B_4 [T@Ref]')
% grid on
% 
% 
% % plot deca
% figure(55)
% plot(traj0(2,:),B_deca,'-r');hold on
% plot(s_alba,b5_alba,'-b')
% xlabel('S (mm)')
% ylabel('Deca B_5 [T@Ref]')
% grid on


%%

% plot grad
figure(2)
plot(traj0(2,:),pn(:,4),'-b');hold off
xlabel('S (mm)')
ylabel('Gradient')
grid on
fprintf('   Sum quad= %8.2e     int %8.2e\n',ps(4), pint(4))



% plot sext
figure(3)
plot(traj0(2,:),pn(:,3),'-b');hold off
xlabel('S (mm)')
ylabel('Sext')
grid on
fprintf('   Sum sext= %8.2e     int %8.2e \n',ps(3),pint(3))

% plot oct
figure(4)
plot(traj0(2,:),pn(:,2),'-b');hold off
xlabel('S (mm)')
ylabel('Oct')
grid on
fprintf('   Sum oct = %8.2e     int %8.2e \n',ps(2),pint(2))

% plot deca
figure(5)
plot(traj0(2,:),pn(:,1),'-b');hold off
xlabel('S (mm)')
ylabel('Deca')
grid on
fprintf('   Sum deca= %8.2e     int %8.2e \n',ps(1),pint(1))


% Output for BETA
% spliting in two parts at 100 mm : inner and outer
lb=0.100;
ne=floor(lb/ds0);
psin=sum(pn(1:ne,:))*ds0/lb; % inner part
psex=ps-psin;                % residual outer part
rx=1.8e-2;                    %radius to get dB/B

fprintf(' \n')
fprintf('  Integrated in: \n')
fprintf('    %8i',[4:-1:0]);fprintf('\n')
fprintf('    %8.2e',psin);fprintf('\n')
fprintf('    %8.2e',psin(1)*(rx)^4,psin(2)*(rx)^3,...
                    psin(3)*(rx)^2);fprintf('\n')
fprintf('  Integrated ex: \n')
fprintf('    %8i',[4:-1:0]);fprintf('\n')
fprintf('    %8.2e',psex);fprintf('\n')
fprintf('    %8.2e',psex(1)*(rx)^4,psex(2)*(rx)^3,...
                    psex(3)*(rx)^2);fprintf('\n')
                
fprintf('\n')
fprintf('Output for BETA : 3 core lens + 2 edge lens \n')
fprintf('\n')

% ps * 2*Lmag (full dipole) / radius / nlens

fprintf('   SXDII   LD  %10.5E 0.6000E+01 \n',psin(3)*2*Lmag/0.352/3)
fprintf('   SXDIE   LD  %10.5E 0.6000E+01 \n',psex(3)*2*Lmag/0.352/2)
fprintf('   OCDII   LD  %10.5E 0.8000E+01 \n',psin(2)*2*Lmag/0.352/3)
fprintf('   OCDIE   LD  %10.5E 0.8000E+01 \n',psex(2)*2*Lmag/0.352/2)
fprintf('   DCDII   LD  %10.5E 0.1000E+02 \n',psin(1)*2*Lmag/0.352/3)
fprintf('   DCDIE   LD  %10.5E 0.1000E+02 \n',psex(1)*2*Lmag/0.352/2)


%  pnvperfect=pn;
%  traj0perfect=traj0;
%  save dipole_perfect.mat traj0perfect pnvperfect

%%

%ps(:,4)*2*Lmag/0.352
%%
% rx=2e-2; 
% 
% DBB_quad = (ps(:,4))*(rx);
% DBB_sext = (ps(:,3))*(rx)^2;
% DBB_oct = (ps(:,2))*(rx)^3;
% DBB_deca = (ps(:,1))*(rx)^4;
% 
% Bint_quad = DBB_quad*bzn(13,1)*2*Lmag;
% Bint_sext = DBB_sext*bzn(13,1)*2*Lmag;
% Bint_oct = DBB_oct*bzn(13,1)*2*Lmag;
% Bint_deca = DBB_deca*bzn(13,1)*2*Lmag;
% 
% fprintf('   Int B2 quad [T m]  %10.5E \n',Bint_quad)
% fprintf('   Int B3 sext [T m]%10.5E  \n',Bint_sext)
% fprintf('   Int B4 octo [T m] %10.5E  \n',Bint_oct)
% fprintf('   Int B5 deca [T m] %10.5E  \n',Bint_deca)

%% focal length

foclengthinv = ps(:,4)*2*Lmag/0.352;

%%

return
 
figure(100); 
plot(xx,polyval(ps,xx*1e-3),'ob');hold on 
plot(xx,polyval(pint,xx*1e-3),'or');hold off 
grid on




