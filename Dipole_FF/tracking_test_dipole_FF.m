
clc; close all; clear all;
%% Initial lattices

thomx_ring=thomx_DIP;
thomx_ring_anglZERO=thomx_DIP_Angle0;


%Z0=[0.005 0 0.005 0 0 0]'*(1:15);
% Z0=[0.01 0 0.01 0 0 0]'*ones(1,15);
Z0=[0.01 0 0.00 0 0 0]'*ones(1,15);

%%

% Z0=zeros(6,15);
% Z0(1,:) = 0.001;

% thomx_ring=thomx;
% thomx_ring_anglZERO=thomx_anglZERO;

[OUT_thomx,lost_thomx]=ringpass(thomx_ring,Z0,1000); %(X, PX, Y, PY, DP, CT2 ) 
[OUT_thomx_anglZERO,lost_anglZERO] =ringpass(thomx_ring_anglZERO,Z0,1000);

 X = OUT_thomx(1,:);
 PX= OUT_thomx(2,:);
 Y= OUT_thomx(3,:);
 PT= OUT_thomx(4,:);
 DP= OUT_thomx(5,:);
 CT= OUT_thomx(6,:);
 
 X_anglZERO = OUT_thomx_anglZERO(1,:);
 PX_anglZERO= OUT_thomx_anglZERO(2,:);
 Y_anglZERO= OUT_thomx_anglZERO(3,:);
 PT_anglZERO= OUT_thomx_anglZERO(4,:);
 DP_anglZERO= OUT_thomx_anglZERO(5,:);
 CT_anglZERO= OUT_thomx_anglZERO(6,:);

figure
plot(X,PX,'b.','DisplayName', 'entrance/Exit angle = 1.9 deg')
hold on
plot(X_anglZERO,PX_anglZERO,'r.','DisplayName', 'entrance/Exit angle = 0')
hold off
legend('show','Location','NorthEast');
%print('tracking','-dpng','-r300')

%% putting FringeBendEntrance/FringeBendExit to 1 
%       FringeBendEntrance/FringeBendExit = 0,1,2,3
%       Version 0 no dipole fringe fields
%       Version 1 legacy version Brown First Order (K. Brown. A First and Second Order 
%                  Matrix Theory for the Design of Beam Transport Systems and Charged 
%                  Particle Spectrometers. Internal report, SLAC-75, 1982)
%       Version 2 SOLEIL close to second order of Brown (J. Bengtsson and M. Meddahi. 
%                 Modeling of Beam Dynamics and Comparison with Measurements for 
%                 the Advanced Light Source. London, UK, 1994.)
%       Version 3 THOMX (Dipole Fringe Field Effects in the ThomX Ring, J. Zhang and 
%                 A. Loulergue, Proceedings of IPAC2013, Shanghai, China)



ring_symplectic =atsetfieldvalues(thomx_ring,findcells(thomx_ring,'PassMethod','BndMPoleSymplectic4Pass'),...
    'FringeBendEntrance',1 );
ring_symplectic =atsetfieldvalues(ring_symplectic,findcells(ring_symplectic,'PassMethod','BndMPoleSymplectic4Pass'),...
    'FringeBendExit',1);


ring_symplectic_anglZERO =atsetfieldvalues(thomx_ring_anglZERO,findcells(thomx_ring_anglZERO,'PassMethod','BndMPoleSymplectic4Pass'),...
    'FringeBendEntrance',1 );
ring_symplectic_anglZERO =atsetfieldvalues(ring_symplectic_anglZERO,findcells(ring_symplectic_anglZERO,'PassMethod','BndMPoleSymplectic4Pass'),...
    'FringeBendExit',1 );



OUT_thomx=ringpass(ring_symplectic,Z0,1000); %(X, PX, Y, PY, DP, CT2 ) 
OUT_thomx_anglZERO=ringpass(ring_symplectic_anglZERO,Z0,1000);

 X = OUT_thomx(1,:);
 PX= OUT_thomx(2,:);
 Y= OUT_thomx(3,:);
 PT= OUT_thomx(4,:);
 DP= OUT_thomx(5,:);
 CT= OUT_thomx(6,:);
 
 X_anglZERO = OUT_thomx_anglZERO(1,:);
 PX_anglZERO= OUT_thomx_anglZERO(2,:);
 Y_anglZERO= OUT_thomx_anglZERO(3,:);
 PT_anglZERO= OUT_thomx_anglZERO(4,:);
 DP_anglZERO= OUT_thomx_anglZERO(5,:);
 CT_anglZERO= OUT_thomx_anglZERO(6,:);

figure
plot(X,PX,'b.','DisplayName', 'entrance/Exit angle = 1.9 deg')
hold on
plot(X_anglZERO,PX_anglZERO,'r.','DisplayName', 'entrance/Exit angle = 0')
hold off
legend('show','Location','NorthEast');
%print('tracking','-dpng','-r300')

%% As you wrote in the email, I am getting NaN if x=y=0.3 which indicates perhpas that particles are lost ???

% Check out this :
% rin=[0.03; 0.00; 0.03; 0; 0; 0];
% [X,lost]=ringpass(RING,rin,1000,'reuse'); 
% figure(50)
% plot(X(1,:),X(2,:),'.b')

ring= thomx_DIP(); 

rin0=[0.005; 0.00; 0.005; 0; 0; 0];
rin=[0.03; 0.00; 0.03; 0; 0; 0];
rin2=[0.01; 0.00; 0.01; 0; 0; 0];
rin3=[0.013; 0.00; 0.013; 0; 0; 0];

[X0,lost0]=ringpass(ring,rin0,1000);
[X,lost]=ringpass(ring,rin,1000); 
[X2,lost2]=ringpass(ring,rin2,1000);
[X3,lost3]=ringpass(ring,rin3,1000);
figure(50)
plot(X(1,:),X(2,:),'.b')

figure(51)
plot(X2(1,:),X2(2,:),'.m')
hold on
plot(X0(1,:),X0(2,:),'.b')
plot(X3(1,:),X3(2,:),'*k')


