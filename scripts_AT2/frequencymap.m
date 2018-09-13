% frequency maps

ring= THOMX();

X=linspace(-0.02,0.02,101);
Y=linspace(-0.0,0.02,51);
nbturn=100;


type=1;

plotnonlinmap(ring,X,Y,nbturn,type,0,'thomxtestfma');


type=2;
plotnonlinmap(ring,X,Y,nbturn,type,0,'thomxtestfma');