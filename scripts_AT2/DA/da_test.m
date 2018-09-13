ring= thomx(); % or THERING

ring0=atfitchrom(ring,[0 0],'SX1\w*','SX3\w*');
ring0=atfitchrom(ring0,[0 0],'SX1\w*','SX3\w*');
ring0=atfitchrom(ring0,[0 0],'SX1\w*','SX3\w*');

%% compute DA
[XX,ZZ]   = atdynapmod(ring ,1000,0); 
[XX0,ZZ0] = atdynapmod(ring0,1000,0); 

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