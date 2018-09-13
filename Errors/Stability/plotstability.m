
%%

stabres = load('stabilityReport.mat');

q= cell2mat(stabres.p.fraction_unstability);
qq = mat2str(q);
frac_unstab = str2num(qq);
%%
figure(1)
set(gca,'FontSize',14)
plot(stabres.p.quaderr, frac_unstab,'.-r','LineWidth',1.3, 'Markersize',15)

xlabel('Quad errors [k0+k0*quaderr*randn]')
ylabel('Fraction of unstable machines');

title(['Number of random machine simulated per QUAD error is: ' num2str(stabres.p.nr_machines) ] );
set(gcf, 'color', 'w');
%export_fig(['stability_QUADerr_' num2str(stabres.p.nr_machines) 'randMach_' num2str(abs(max(stabres.p.quaderr))) 'K.pdf'])

%%
figure(2)
set(gca,'FontSize',14)
stem(stabres.p.quaderr, frac_unstab,'LineStyle','-.','LineWidth',1.3,...
     'MarkerFaceColor','red',...
     'MarkerEdgeColor','green', 'Markersize',7)

xlabel('Quad errors [k0+k0*quaderr*randn]')
ylabel('Fraction of unstable machines');

title(['Number of random machine simulated per QUAD error is: ' num2str(stabres.p.nr_machines) ] );
set(gcf, 'color', 'w');
%export_fig(['stability2_QUADerr_' num2str(stabres.p.nr_machines) 'randMach_' num2str(abs(max(stabres.p.quaderr))) 'K.pdf'])