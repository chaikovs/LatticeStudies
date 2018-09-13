close all; clear all

r_ml = thomxML();

%%

r_optic=matchthomx_diplength(r_ml);


%%
figure; c=atplot(r_ml);
figure; c=atplot(r_optic);   % after match
%%

k_original = getkval(r_ml);
k_match = getkval(r_optic);
%%

% plot 
% figure; c=atplot(r_original);  c.comment.String{3}=[' alpha= ' num2str(mcf(r))];% before match
% figure; c=atplot(r_ml);  c.comment.String{3}=[' alpha= ' num2str(mcf(r_optic))]; % after match

figure; c=atplot(r_ml);  %print('beta_atplot_ML.png', '-dpng', '-r300');% before match
figure; c=atplot(r_optic); %print('beta_atplot_Match.png', '-dpng', '-r300');  % after match
%%

figure
subplot 211
set(gca,'FontSize',18)
plot(k_match(1,:),'ro-','DisplayName', 'After matching')
hold on
plot(k_original(1,:),'bo-','DisplayName', 'Before matching')
hold off
legend('show','Location','NorthWest');
xlabel('Family number');
ylabel('K value [1/m^2]');
subplot 212
set(gca,'FontSize',18)
bar(k_match(1,:)./k_original(1,:))
xlabel('Family number');
ylabel('K_{match}/K_{orig}');
print('K_matching.png', '-dpng', '-r300');
%%
%indBPM=find(atgetcells(r_optic,'FamName','BPMx'))';
ind=1:length(r_optic);

[l,t,c]=atlinopt(r_optic,0,ind);
[le,te,ce]=atlinopt(r_ml,0,ind);
bx0=arrayfun(@(a)a.beta(1),le);
bxe=arrayfun(@(a)a.beta(1),l);


figure;
plot(bx0); 
hold all; 
plot(bxe); 
hold off
legend('before match','after match')
title([' Before matching ' num2str(te) '; After matching ' num2str(t)])

%%

% qp1ind     = findcells(ring,'FamName','QP1');
% qp2ind     = findcells(ring,'FamName','QP2');
% qp3ind     = findcells(ring,'FamName','QP3');
% qp4ind     = findcells(ring,'FamName','QP4');
% qp31ind     = findcells(ring,'FamName','QP31');
% qp41ind     = findcells(ring,'FamName','QP41');

atind = atindex(r_ml);


r_ml_new = atsetfieldvalues(r_ml,atind.QP1','PolynomB',{1,2},k_match(:,1)); %atsetfieldvalues
r_ml_new = atsetfieldvalues(r_ml_new,atind.QP2,'PolynomB',{1,2},k_match(:,2));
r_ml_new = atsetfieldvalues(r_ml_new,atind.QP3 ,'PolynomB',{1,2},k_match(:,3));
r_ml_new = atsetfieldvalues(r_ml_new,atind.QP4,'PolynomB',{1,2},k_match(:,4));
r_ml_new = atsetfieldvalues(r_ml_new,atind.QP31,'PolynomB',{1,2},k_match(:,5));
r_ml_new = atsetfieldvalues(r_ml_new,atind.QP41,'PolynomB',{1,2},k_match(:,6));

[~,tune,~]=atlinopt(r_ml,0,1)
[~,tune_new,~]=atlinopt(r_ml_new,0,1)

%%
