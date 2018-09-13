
clear all; close all; clc

%%
ring_DIPmagnL= thomx_DIP();  
ring_DIPmagnL_Angle0 = thomx_DIP_Angle0();
ring_DIP_QUADmagnL= thomx_QAUD();

beta_DA_DIPmagnL = load('DA_DIPmagnL');
beta_DA_DIPmagnL_Angle0 = load('DA_DIP_Angle0');
beta_DA_DIP_QUADmagnL = load('DA_DIP_QUADMagnL');

ring_DIPmagnL_linear=atsetfieldvalues(ring_DIPmagnL,findcells(ring_DIPmagnL,'PassMethod','BndMPoleSymplectic4Pass'),...
    'PassMethod','BendLinearPass');

ring_DIPmagnL_linear_chro0=atfitchrom(ring_DIPmagnL_linear,[0 0],'SX1\w*','SX2\w*');
ring_DIPmagnL_linear_chro0=atfitchrom(ring_DIPmagnL_linear_chro0,[0 0],'SX1\w*','SX2\w*');
ring_DIPmagnL_linear_chro0=atfitchrom(ring_DIPmagnL_linear_chro0,[0 0],'SX1\w*','SX2\w*');


ring_DIPmagnL_chro0=atfitchrom(ring_DIPmagnL,[0 0],'SX1\w*','SX2\w*');
ring_DIPmagnL_chro0=atfitchrom(ring_DIPmagnL_chro0,[0 0],'SX1\w*','SX2\w*');
ring_DIPmagnL_chro0=atfitchrom(ring_DIPmagnL_chro0,[0 0],'SX1\w*','SX2\w*');


%% compute DA


[XX_DIPmagnL,ZZ_DIPmagnL]   = atdynap(ring_DIPmagnL, 1000,0); 
[XX_DIPmagnL_linear,ZZ_DIPmagnL_linear]   = atdynap(ring_DIPmagnL_linear, 1000,0); 

[XX_DIPmagnL_chro0, ZZ_DIPmagnL_chro0]   = atdynap(ring_DIPmagnL_chro0, 1000,0);
[XX_DIPmagnL_linear_chro0,ZZ_DIPmagnL_linear_chro0]   = atdynap(ring_DIPmagnL_linear_chro0, 1000,0); 

[XX_DIPmagnL_Angle0,ZZ_DIPmagnL_Angle0]   = atdynap(ring_DIPmagnL_Angle0, 1000,0); 

[XX_QUADmagnL,ZZ_QUADmagnL] = atdynap(ring_DIP_QUADmagnL,1000,0); 


%% plot DA

figure;
set(gca,'FontSize',20)
plot(XX_DIPmagnL,ZZ_DIPmagnL,'r.-'); 
hold on; 
plot(beta_DA_DIPmagnL(:,1),beta_DA_DIPmagnL(:,2),'g-')
plot(XX_DIPmagnL_chro0,ZZ_DIPmagnL_chro0,'b.-'); 
xlabel('x')
ylabel('y')

[l_DIPmagnL,t_DIPmagnL,c_DIPmagnL]=atlinopt(ring_DIPmagnL,0,1);
[l_DIPmagnL_chro0,t_DIPmagnL_chro0,c_DIPmagnL_chro0]=atlinopt(ring_DIPmagnL_chro0,0,1);


legend(['Lattice DIP magnL ' num2str(c_DIPmagnL)], 'Lattice DIP magnL BETA code', ['Lattice DIP magnL Chro 0 ' num2str(c_DIPmagnL_chro0)])
print('DA_thomx2','-dpng','-r300')

%%

figure;
set(gca,'FontSize',20)
plot(XX_DIPmagnL,ZZ_DIPmagnL,'r.-'); 
hold on; 
plot(beta_DA_DIPmagnL(:,1),beta_DA_DIPmagnL(:,2),'g--')
plot(XX_DIPmagnL_linear,ZZ_DIPmagnL_linear,'b.-'); 
plot(XX_DIPmagnL_linear_chro0,ZZ_DIPmagnL_linear_chro0,'m.-'); 
xlabel('x')
ylabel('y')

[l_DIPmagnL,t_DIPmagnL,c_DIPmagnL]=atlinopt(ring_DIPmagnL,0,1);
[l_DIPmagnL_chro0,t_DIPmagnL_chro0,c_DIPmagnL_chro0]=atlinopt(ring_DIPmagnL_chro0,0,1);
[l_DIPmagnL_linear_chro0,t_DIPmagnL_linear_chro0,c_DIPmagnL_linear_chro0]=atlinopt(ring_DIPmagnL_linear_chro0,0,1);


legend(['Lattice DIP magnL ' num2str(c_DIPmagnL)], 'Lattice DIP magnL BETA code', ['Lattice DIP magnL BendLinearPass ' num2str(c_DIPmagnL_chro0)], ['Lattice DIP magnL BendLinearPass Chro 0 ' num2str(c_DIPmagnL_linear_chro0)])
print('DA_thomx','-dpng','-r300')

%%

figure;
set(gca,'FontSize',20)
plot(XX_DIPmagnL,ZZ_DIPmagnL,'b.-'); 
hold on; 
plot(beta_DA_DIPmagnL(:,1),beta_DA_DIPmagnL(:,2),'b--')
plot(XX_QUADmagnL,ZZ_QUADmagnL,'r.-'); 
plot(beta_DA_DIP_QUADmagnL(:,1),beta_DA_DIP_QUADmagnL(:,2),'r--')
xlabel('x')
ylabel('y')

[l_DIPmagnL,t_DIPmagnL,c_DIPmagnL]=atlinopt(ring_DIPmagnL,0,1);
[l_DIP_QUADmagnL,t_DIP_QUADmagnL,c_DIP_QUADmagnL]=atlinopt(ring_DIP_QUADmagnL,0,1);


legend(['Lattice DIP magnL ' num2str(c_DIPmagnL)], 'Lattice DIP magnL BETA code', ['Lattice DIP QUAD magnL ' num2str(c_DIP_QUADmagnL)], 'Lattice DIP QUAD magnL BETA code')


%% plot DA

figure;
set(gca,'FontSize',20)
plot(XX_DIPmagnL,ZZ_DIPmagnL,'b.-'); 
hold on; 
plot(beta_DA_DIPmagnL(:,1),beta_DA_DIPmagnL(:,2),'b--')
plot(XX_DIPmagnL_Angle0,ZZ_DIPmagnL_Angle0,'r.-'); 
plot(beta_DA_DIPmagnL_Angle0(:,1),beta_DA_DIPmagnL_Angle0(:,2),'r--')
xlabel('x')
ylabel('y')

[l_DIPmagnL,t_DIPmagnL,c_DIPmagnL]=atlinopt(ring_DIPmagnL,0,1);
[l_DIPmagnL_Angle0,t_DIPmagnL_Angle0,c_DIPmagnL_Angle0]=atlinopt(ring_DIPmagnL_Angle0,0,1);


legend(['Lattice DIP magnL ' num2str(c_DIPmagnL)], 'Lattice DIP magnL BETA code', ['Lattice DIP Angle 0 magnL ' num2str(c_DIPmagnL_Angle0)], 'Lattice DIP Angle 0 BETA code')
print('DA_thomx_Angle0','-dpng','-r300')


