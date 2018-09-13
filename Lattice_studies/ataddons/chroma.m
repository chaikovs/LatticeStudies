%
clear
DP=-0.04;
%
dir ='/home/sources/physmach/loulergue/work/matlab/Simulation/beta_structure/';
%input_file =[dir 'Upgrade/test-offmomentum/6BA-DLS-3Q3-fit.str'];
input_file =[dir 'Upgrade/6BA-7BA-AB/Full-ring-2/6BA-DLS-scaled-2-sx2-fit1-optim17.str'];
%ring0=atreadbeta_alex(input_file);

%load('Lattice-test/lat_7BA_ESRF_3Q_atmatch40');
load('Lattice-save/SOLEIL-U-v4');
ring0=RING;

%
%ring0=fitbend(ring0,0.3993,[0.72 1.02],'6BA');
ring0=scalesext(ring0,'SXD1E',1.);
ring0=fitchrom_alex(ring0,[0.2 0],'SXD2E' ,'SXF1E' );
%
[px ,pz]=get_nudp(DP,ring0);
fprintf('px =  %8.2f   %8.2f   %8.2f   %8.2f  \n',px)
fprintf('pz =  %8.2f   %8.2f   %8.2f   %8.2f  \n',pz)
%
el1=1;
el2=length(ring0)+1;
totlength=findspos(ring0(el1:el2-1),length(ring0(el1:el2-1))+1);
ring=atslice(ring0(el1:el2-1),totlength*10);
el3=length(ring)+1;

spos  = findspos(ring,1:el3);
[lindata,tunes,chrom]=atlinopt(ring,0,1:el3); 
disp  = cat(2,lindata.Dispersion);
phase0= cat(1,lindata.mu)';
beta0 = cat(1,lindata.beta)';

orbit = findorbit4(ring,DP,1:el3);
[nonlindata] = twissring(ring,DP,1:el3, 'chrom');
phase = cat(1,nonlindata.mu)';
beta  = cat(1,nonlindata.beta)';
%
figure(3)
set(gcf,'color','w')
h1=subplot(4,1,[1 3]);
plot(spos,disp(1,:)*DP*1e3,'-k'); hold on
plot(spos,orbit(1,:)*1e3,'-r');
plot(spos,(orbit(1,:)-disp(1,:)*DP)*1e3,'-m'); hold off
xlim([0 spos(end)]);
ylabel('X (mm)')
legend('Linear','Non-linear','Diff')
h2 = subplot(4,1,4);
atplotsyn(h2,ring0);
xlabel('S (m)')
linkaxes([h1 h2],'x')
set(h1,'xtick',[]);
set(h2,'ytick',[]);

return
figure(31)
set(gcf,'color','w')
h1=subplot(4,1,[1 3]);
plot(spos,beta0(1,:),'-k'); hold on
plot(spos,beta(1,:),'-r');
plot(spos,(beta(1,:)-beta0(1,:))./beta0(1,:)*100,'-m'); hold off
xlim([0 spos(end)]);
ylabel('X (mm)')
legend('Linear','Non-linear','Diff relative')
h2 = subplot(4,1,4);
atplotsyn(h2,ring0);
xlabel('S (m)')
linkaxes([h1 h2],'x')
set(h1,'xtick',[]);
set(h2,'ytick',[]);

figure(32)
set(gcf,'color','w')
h1 = subplot(4,1,[1 3]);
plot(spos,beta0(2,:),'-k'); hold on
plot(spos,beta(2,:),'-r');
plot(spos,(beta(2,:)-beta0(2,:))./beta0(2,:)*100,'-m'); hold off
xlim([0 spos(end)]);
ylabel('X (mm)')
legend('Linear','Non-linear','Diff relative')
h2 = subplot(4,1,4);
atplotsyn(h2,ring0);
xlabel('S (m)')
linkaxes([h1 h2],'x')
set(h1,'xtick',[]);
set(h2,'ytick',[]);

%
figure(4)
set(gcf,'color','w')
h1 = subplot(4,1,[1 3]);
plot(spos,(phase(1,:)-phase0(1,:))/2/pi,'-r'); hold on
plot(spos,(phase(2,:)-phase0(2,:))/2/pi,'-b'); hold off
xlim([0 spos(end)]);
ylabel('Phase')
h2 = subplot(4,1,4);
atplotsyn(h2,ring0);
xlabel('S (m)')
linkaxes([h1 h2],'x')
set(h1,'xtick',[]);
set(h2,'ytick',[]);



