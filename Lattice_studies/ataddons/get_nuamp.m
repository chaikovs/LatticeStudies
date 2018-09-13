function [Q] = get_nuamp(A,RING,order,dp,c,nfig)
% get nuamp versus radial amplitude at fixed dp
[~, tunes] = twissring(RING, 0, 1:length(RING)+1);
period   = 1;
if nargin>=5;period=c;end

intgtunex=floor(tunes(1));
intgtunez=floor(tunes(2));
fractunex=tunes(1)-intgtunex;
fractunez=tunes(2)-intgtunez;
%
ord=2;if nargin > 2 ; ord=order ;end
n=7;
x=(-A(1) : A(1)/n : A(1));
z=(-A(2) : A(2)/n : A(2));
%
[nuxx,nuzx]=nuampl_alex(RING,x,1,dp,0);
if (fractunex>0.5); nuxx=1-nuxx;end
if (fractunez>0.5); nuzx=1-nuzx;end
nuxx=intgtunex+nuxx;
nuzx=intgtunez+nuzx;
px = polyfit_nowarning(x,nuxx,ord);Q(1)=px(ord-1)/1;
pz = polyfit_nowarning(x,nuzx,ord);Q(2)=pz(ord-1)/1;
fnuxx = polyval(px,x); 
fnuzx = polyval(pz,x); 
% 
%ord=2;% enought for vertical
[nuxz,nuzz]=nuampl_alex(RING,z,3,dp,0);
if (fractunex>0.5); nuxz=1-nuxz;end
if (fractunez>0.5); nuzz=1-nuzz;end
nuxz=intgtunex+nuxz;
nuzz=intgtunez+nuzz;
px = polyfit(z,nuxz,ord);Q(3)=px(ord-1)/1;
pz = polyfit(z,nuzz,ord);Q(4)=pz(ord-1)/1;
fnuxz = polyval(px,z); 
fnuzz = polyval(pz,z); 
%

f(1)=15;f(2)=16;f(3)=17;
if(nargin==6);f=nfig;end

% 
figure(f(1))
set(gcf,'color','w');
plot(x*1e3,(nuxx-nuxx(n+1))*period,'or','LineWidth',2);hold on
plot(x*1e3,(nuzx-nuzx(n+1))*period,'ob','LineWidth',2);
plot(x*1e3,(fnuxx-nuxx(n+1))*period,'-r','LineWidth',2);
plot(x*1e3,(fnuzx-nuzx(n+1))*period,'-b','LineWidth',2);hold off
grid on
set(gca,'fontsize',20);
xlabel('x [mm]');               
ylabel('\delta\nu');
legend('X','Z')
title('Tunes vs x')
print('thomx_lattice_nuamp_nux_WP0_multip.png','-dpng','-r300')

figure(f(2))
set(gcf,'color','w');
plot(z*1e3,(nuxz-nuxz(n+1))*period,'or','LineWidth',2);hold on
plot(z*1e3,(nuzz-nuzz(n+1))*period,'ob','LineWidth',2);
plot(z*1e3,(fnuxz-nuxz(n+1))*period,'-r','LineWidth',2);
plot(z*1e3,(fnuzz-nuzz(n+1))*period,'-b','LineWidth',2);hold off
grid on
set(gca,'fontsize',20);
xlabel('z [mm]');               
ylabel('\delta\nu');
legend('X','Z')
title('Tunes vs z')
print('thomx_lattice_nuamp_nuz_WP0_multip.png','-dpng','-r300')

cnx=0;cnz=0;
intgtunex=floor(tunes(1)*period);
intgtunez=floor(tunes(2)*period);
fractunex=tunes(1)*period-intgtunex;
fractunez=tunes(2)*period-intgtunez;
if (fractunex>0.5); cnx=0.5; end
if (fractunez>0.5); cnz=0.5; end

figure(f(3))
title('Tunes vs \delta')
set(gcf,'color','w');
plot(nuxx*period,nuzx*period,'or','LineWidth',2);hold on
plot(nuxz*period,nuzz*period,'ob','LineWidth',2);
plot(nuxx(n+1)*period,nuzx(n+1)*period,'ok','MarkerSize',5,'MarkerFaceColor','k');hold off
xlim([cnx cnx+0.5]+intgtunex)
ylim([cnz cnz+0.5]+intgtunez)
set(gca,'fontsize',20);
xlabel('\nu_x');
ylabel('\nu_z');
legend('X','Z')
grid on
print('thomx_lattice_nuamp_nuxnuz_WP0_multip.png','-dpng','-r300')
end

