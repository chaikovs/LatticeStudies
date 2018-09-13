
%%

ring= thomx_FF(); 


%%

indx=1:length(ring);    
T=twissring(ring,0,indx);
beta=cat(1,T.beta);
   
%%
rx_bpipe = 20e-3;
rz_bpipe = 14e-3;

bxinj = max(beta(1,1));
bzinj = max(beta(1,2));
bxmax = max(beta(:,1));
bzmax = max(beta(:,2));

rx_bpipe_scaled = rx_bpipe / sqrt(bxmax/bxinj)
rz_bpipe_scaled = rz_bpipe / sqrt(bzmax/bzinj)

a=rx_bpipe_scaled; % horizontal radius
b=rz_bpipe_scaled; % vertical radius
x0=0; % x0,y0 ellipse centre coordinates
y0=0;
t=-pi:0.01:pi;
x=x0+a*cos(t);
y=y0+b*sin(t);

a_bpipe=rx_bpipe; % horizontal radius
b_bpipe=rz_bpipe; % vertical radius
x0_bpipe=0; % x0,y0 ellipse centre coordinates
y0_bpipe=0;
t_bpipe=-pi:0.01:pi;
x_bpipe=x0_bpipe+a_bpipe*cos(t_bpipe);
y_bpipe=y0_bpipe+b_bpipe*sin(t_bpipe);

%%

[XX,ZZ]   = atdynap(ring, 1000,0.0,0.02);
[XX_off_minus2,ZZ_off_minus2]   = atdynap(ring, 1000,-0.02, 0.02);
[XX_off_plus2,ZZ_off_plus2]   = atdynap(ring, 1000, 0.02, 0.02);
[XX_off_minus15,ZZ_off_minus15]   = atdynap(ring, 1000,-0.015, 0.02);
[XX_off_plus15,ZZ_off_plus15]   = atdynap(ring, 1000, 0.015, 0.02);

%%

figure('units','normalized','position',[0.3 0.3 0.45 0.35])
plot(x,y,'k-','DisplayName', 'Scaled vacuum chamber')
hold on; 
plot(x_bpipe,y_bpipe,'k--','DisplayName', 'Vacuum chamber')
plot(XX,ZZ,'go-','MarkerSize',7,'DisplayName', 'DA ON-momentum'); %,'HandleVisibility','off' 
return
plot(XX_off_minus2,ZZ_off_minus2,'r*-','MarkerSize',10,'DisplayName', 'DA OFF-momentum -2%');
plot(XX_off_plus2,ZZ_off_plus2,'rd-','MarkerSize',10,'DisplayName', 'DA OFF-momentum +2%');
plot(XX_off_plus15,ZZ_off_plus15,'bd-','MarkerSize',10,'DisplayName', 'DA OFF-momentum +1.5%');
plot(XX_off_minus15,ZZ_off_minus15,'b*-','MarkerSize',10,'DisplayName', 'DA OFF-momentum -1.5%');
xlabel('x [m]')
ylabel('z [m]')
set(gca,'FontSize',18)
set(gcf,'color','w')
u = legend('show','Location','NorthEast');
set(u,'FontSize',14)
xlim([-0.04 0.04])
ylim([0 0.025])
print('DA_OFFmomentum','-dpng','-r300')

%%


