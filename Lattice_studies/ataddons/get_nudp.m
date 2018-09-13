function [px ,pz]=get_nudp(d, RING, c)
% get nudp
global THERING
%
if nargin>=2;
    ring0=RING;
else
    ring0=THERING;
end
%
ncell=1;
if nargin==3;ncell=c;end
%
dp=(-d : d/10 : d);
i=0;
for dd=dp
    i=i+1;
    %[TD, tune] = twissring(THERING,dd,1);
    [~, sum.tunes] = twissring(ring0, dd, 1:length(ring0)+1);
    nux(i)=sum.tunes(1);nuz(i)=sum.tunes(2);
end

nux=nux*ncell;
nuz=nuz*ncell;

% Get polynome
px = polyfit(dp,nux,4);
fnux = polyval(px,dp); 
pz = polyfit(dp,nuz,4);
fnuz = polyval(pz,dp); 

n0=11;
% Plot
figure(1)
set(gcf,'color','w');
plot(dp*1e2,nux-nux(n0),'or','LineWidth',2);hold on
plot(dp*1e2,nuz-nuz(n0),'ob','LineWidth',2);
plot(dp*1e2,fnux-nux(n0),'-r','LineWidth',2);
plot(dp*1e2,fnuz-nuz(n0),'-b','LineWidth',2);hold off
grid on
set(gca,'fontsize',20);
xlabel('\deltap [%]');                 % Add labels
ylabel('\delta\nu');
legend('X','Z')
title('Tunes vs \delta')
%xlim([-3 3])
%ylim([-0.06 0.06])
print('thomx_lattice_dnudp_WP0_multip.png','-dpng','-r300')

figure(11)
set(gcf,'color','w');
plot(dp*1e2,nux,'or','LineWidth',2);hold on
plot(dp*1e2,nuz,'ob','LineWidth',2);
plot(dp*1e2,fnux,'-r','LineWidth',2);
plot(dp*1e2,fnuz,'-b','LineWidth',2);hold off
grid on
set(gca,'fontsize',20);
xlabel('\deltap [%]');                 % Add labels
ylabel('\delta\nu');
legend('X','Z')
title('Tunes vs \delta')
%xlim([-3 3])
%ylim([-0.06 0.06])
print('thomx_lattice_dnudp_WP0_tune_multip.png','-dpng','-r300')

% Plot
cnx=0;cnz=0;
if (nux(n0)-floor(nux(n0))>0.5) ; cnx=0.5; end
if (nuz(n0)-floor(nuz(n0))>0.5) ; cnz=0.5; end

figure(2)
title('Tunes vs \delta')
set(gcf,'color','w');

plot(nux,nuz,'ok','LineWidth',2);hold on
%plot(fnux,fnuz,'-k','LineWidth',2);
plot(nux(n0),nuz(n0),'or','MarkerSize',5,'MarkerFaceColor','r');hold off
% xlim([floor(nux(n0))+0.25 floor(nux(n0))+0.5])
% ylim([floor(nuz(n0))+0.75 floor(nuz(n0))+1])
% xlim([floor(nux(n0))+cnx+0.14 floor(nux(n0))+cnx+0.2])
% ylim([floor(nuz(n0))+cnz*1.22 floor(nuz(n0))+cnz+0.2])
%xlim([floor(nux(n0))+cnx floor(nux(n0))+cnx+0.5])
%ylim([floor(nuz(n0))+cnz floor(nuz(n0))+cnz+0.5])
set(gca,'fontsize',20);
xlabel('\nu_x');
ylabel('\nu_z');
grid on
print('thomx_lattice_nuxnuz_WP0_multip.png','-dpng','-r300')


