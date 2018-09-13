close all

r=thomx_DIP_Angle002();

r_optic=matchthomx(r);

% plot 
figure; c=atplot(r);  c.comment.String{3}=[' alpha= ' num2str(mcf(r))];% before match
figure; c=atplot(r_optic);  c.comment.String{3}=[' alpha= ' num2str(mcf(r_optic))]; % after match

%%

r_optic2=matchthomx(r_optic);

% plot 
figure; c=atplot(r_optic);  c.comment.String{3}=[' alpha= ' num2str(mcf(r_optic))];% before match
figure; c=atplot(r_optic2);  c.comment.String{3}=[' alpha= ' num2str(mcf(r_optic2))]; % after match

%%

r_optic3=matchthomx(r_optic2);

% plot 
figure; c=atplot(r_optic2);  c.comment.String{3}=[' alpha= ' num2str(mcf(r_optic2))];% before match
figure; c=atplot(r_optic3);  c.comment.String{3}=[' alpha= ' num2str(mcf(r_optic3))]; % after match


%%

r_optic4=matchthomx(r_optic3);

% plot 
figure; c=atplot(r_optic3);  c.comment.String{3}=[' alpha= ' num2str(mcf(r_optic3))];% before match
figure; c=atplot(r_optic4);  c.comment.String{3}=[' alpha= ' num2str(mcf(r_optic4))]; % after match

%%

r_optic5=matchthomx(r_optic4);

% plot 
figure; c=atplot(r_optic4);  c.comment.String{3}=[' alpha= ' num2str(mcf(r_optic4))];% before match
figure; c=atplot(r_optic5);  c.comment.String{3}=[' alpha= ' num2str(mcf(r_optic5))]; % after match

%%

r_optic6=matchthomx(r_optic5);

% plot 
figure; c=atplot(r_optic5);  c.comment.String{3}=[' alpha= ' num2str(mcf(r_optic4))];% before match
figure; c=atplot(r_optic6);  c.comment.String{3}=[' alpha= ' num2str(mcf(r_optic6))]; % after match

