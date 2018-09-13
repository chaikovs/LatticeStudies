close all

r=thomx02();

r_optic=matchthomx(r);

% plot 
figure; c=atplot(r);  c.comment.String{3}=[' alpha= ' num2str(mcf(r))];% before match
figure; c=atplot(r_optic);  c.comment.String{3}=[' alpha= ' num2str(mcf(r_optic))]; % after match

