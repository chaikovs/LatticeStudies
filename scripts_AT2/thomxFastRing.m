tx=thomxBnd;
[r0 rrad]=atfastring(tx);  
%r0 now contains fast ring with no radiation damping/diffusion
%rrad contains fast ring with radiation effects


%compare dynamic apertures
tic
[x,y]=atdynap(tx,200,0);
toc

tic
[xf,yf]=atdynap(r0,200,0);
toc

figure
plot(x,y,'-r',xf,yf,'-b')

%damping and diffusion

%check symplecticity
m66=findm66(tx);
symperr=m66'*jmat(3)*m66-jmat(3);
symperr


nturns=66e6;
x1=.001;
z1=[x1;0;0;0;0;0];

z0=[0;0;0;0;0;0];

zf0=ringpass(rrad,z0,nturns);
zf1=ringpass(rrad,z1,nturns);
close all
figure
plot(zf0(1,1:1000:end),'-r');
hold on
plot(zf1(1,1:1000:end),'-b');
