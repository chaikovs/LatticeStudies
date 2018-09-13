ring= thomx02(); % or THERING

ring0=atfitchrom(ring,[0 0],'SX1\w*','SX3\w*');
ring0=atfitchrom(ring0,[0 0],'SX1\w*','SX3\w*');
ring0=atfitchrom(ring0,[0 0],'SX1\w*','SX3\w*');

%% compute DA
[XX,ZZ]   = atdynap(ring ,1000,0); 
[XX0,ZZ0] = atdynap(ring0,1000,0); 

%% plot DA
figure;
plot(XX,ZZ); 
hold on; 
plot(XX0,ZZ0);
xlabel('x')
ylabel('y')

[l,t,c]=atlinopt(ring,0,1);
[l0,t0,c0]=atlinopt(ring0,0,1);

legend(['initial ' num2str(c)],['chroma 0,0 ' num2str(c0)])

%% compute lifetime
ring=atsetfieldvalues(ring,findcells(ring,'PassMethod','BendLinearPass'),...
    'PassMethod','BndMPoleSymplectic4Pass');
TL=TouschekPiwinskiLifeTime(...
 ring,...
 0.01,...  % column array (size of r or positions)
 16.7e-3,...                 % scalar
 1:length(ring),...  %(default all elements with length>0 )  column array
 50e-9,... %(default atx modemittance(1))   scalar
 50e-9,... %(default 10 pm)		       scalar
 'trapz',...  % (default quad, may be: 'integral', 'quad', 'trapz', 'elegantLike', 'Approximate')
 0.001,...	% scalar
 0.0075...        % scalar	   
 );
Tl
