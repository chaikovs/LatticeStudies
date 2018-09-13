ring= thomx_DIP(); 
ring_q= thomx_QAUD(); 
ring_AngleZero= thomx_DIP_Angle0(); 

ring_symplectic =atsetfieldvalues(ring,findcells(ring,'PassMethod','BndMPoleSymplectic4Pass'),...
    'FringeBendEntrance',1 );
ring_symplectic =atsetfieldvalues(ring_symplectic,findcells(ring_symplectic,'PassMethod','BndMPoleSymplectic4Pass'),...
    'FringeBendExit',1);

ring0=atfitchrom(ring,[0 0],'SX1\w*','SX2\w*');
ring0=atfitchrom(ring0,[0 0],'SX1\w*','SX2\w*');
ring0=atfitchrom(ring0,[0 0],'SX1\w*','SX2\w*');

ring_q0=atfitchrom(ring_q,[0 0],'SX1\w*','SX2\w*');
ring_q0=atfitchrom(ring_q0,[0 0],'SX1\w*','SX2\w*');
ring_q0=atfitchrom(ring_q0,[0 0],'SX1\w*','SX2\w*');

%% compute DA
[XX,ZZ]   = atdynap(ring, 1000,0); 
[XX0,ZZ0]   = atdynap(ring0, 1000,0); 
 
[XXq,ZZq]   = atdynap(ring_q ,1000,0); 
[XXq0,ZZq0]   = atdynap(ring_q0 ,1000,0);

[XX_ang0,ZZ_ang0] = atdynap(ring_AngleZero,1000,0); 

[XX_symplectic,ZZ_symplectic] = atdynap(ring_symplectic,1000,0); 

%% plot DA
figure;
set(gca,'FontSize',20)
plot(XX,ZZ,'bo-'); 
hold on; 
plot(XX0,ZZ0,'b*-');
plot(XX_symplectic,ZZ_symplectic,'r*-');
plot(XX_ang0,ZZ_ang0,'ro-');
plot(XXq,ZZq,'mo-');
plot(XXq0,ZZq0,'m*-');
xlabel('x')
ylabel('y')

[l,t,c]=atlinopt(ring,0,1);
[l0,t0,c0]=atlinopt(ring0,0,1);
[lq,tq,cq]=atlinopt(ring_q,0,1);
[lq0,tq0,cq0]=atlinopt(ring_q0,0,1);
[l_ang0,t_ang0,c_ang0]=atlinopt(ring_AngleZero,0,1);
[l_symplectic,t_symplectic,c_symplectic]=atlinopt(ring_symplectic,0,1);


legend(['Old lattice ' num2str(c)],['Old lattice Chro 0 ' num2str(c0)],['Old lattice Symplectic ' num2str(c_symplectic)],...
   ['Old lattice Angle 0 ' num2str(c_ang0)], ['New lattice ' num2str(cq)], ['New lattice Chro 0 ' num2str(cq0)])



