function [pa]=plot_momentum_compaction(RING,delta)
%UNTITLED Summary of this function goes here
%   mcf versus delta

dd=-delta:delta/10:delta;
i=0;
for d=dd;
    i =i+1;
    %a(i) = mcf2(RING,d); % second order ?
    a(i) = mcf(RING,d);
end
pa = polyfit(dd,a,3);

figure(1)
set(gcf,'color','w');
set(gca,'fontsize',16);
plot(dd*1e2,a,'ob','LineWidth',2);%hold on
%plot(dd*1e2,fnuz-nuz(n0),'-b','LineWidth',2);hold off
grid on
xlabel('\delta [%]');               
ylabel('MCF');


