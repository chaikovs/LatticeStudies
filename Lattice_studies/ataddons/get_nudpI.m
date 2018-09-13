function [px ,pz]=get_nudpI(d, RING)
% get nudp
global THERING
%
id=findcells(RING,'FamName','SXF1E');
n1=id(1);n2=id(2); 

if nargin==2;
    ring0=RING;
else
    ring0=THERING;
end
%
dp=(-d : d/10 : d);
i=0;
for dd=dp
    i=i+1;
    [TD] = twissring(ring0, dd, 1:length(ring0)+1);
    nu    = cat(1, TD.mu)/2/pi;
    nux(i)=nu(n2,1)-nu(n1,1);
    nuz(i)=nu(n2,2)-nu(n1,2);
end

% Get polynome
px = polyfit(dp,nux,3);
fnux = polyval(px,dp); 
pz = polyfit(dp,nuz,3);
fnuz = polyval(pz,dp); 

n0=11;
% Plot
figure(3)
set(gcf,'color','w');
set(gca,'fontsize',16);
plot(dp*1e2,nux-nux(n0),'or','LineWidth',2);hold on
plot(dp*1e2,nuz-nuz(n0),'ob','LineWidth',2);
plot(dp*1e2,fnux-nux(n0),'-r','LineWidth',2);
plot(dp*1e2,fnuz-nuz(n0),'-b','LineWidth',2);hold off
grid on
xlabel('\delta [%]');                 % Add labels
ylabel('\delta\nu');
legend('X','Z')
title('Tunes vs \delta')
xlim([-6 6])
ylim([-0.06 0.06])


% Plot
figure(4)
title('Tunes vs \delta')
set(gcf,'color','w');
set(gca,'fontsize',16);
plot(nux,nuz,'ok','LineWidth',2);hold on
plot(fnux,fnuz,'-k','LineWidth',2);hold off
grid on


