
%%

ring = ThomX_017_058_r56_02_chro22_multip();

%%

[xmaxlist1,dplist1] = atdynap_om(rerr,(0:2:38)*1e-3,1e-12,(-3:0.2:3)*1e-2,500); 

%%

figure(13);
set(gcf,'color','w')
plot(dplist1*1e2,xmaxlist1*1e3,'-b','LineWidth',2,'DisplayName', 'DA OFF-momentum 3.17/1.64 Chro 0/0');
xlim([-3 3])
ylim([0 max(xmaxlist1)]*1.4e3)
grid on
set(gca,'fontsize',20);
%title('OFF-momentum DA')
xlabel('\deltap [%]');                 % Add labels
ylabel('x [mm]');
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)