
%%

ring= ThomX_017_064_r56_02(); 

%td = load('TuneDiagram');
td1 = load('TuneDiagram_1');
td2 = load('TuneDiagram_2');
td3 = load('TuneDiagram_3');
td4 = load('TuneDiagram_4');
td5 = load('TuneDiagram_5');


[l,t,c] = atlinopt(ring,0,1);

% figure
% set(gcf,'color','w')
% plot(td(:,1), td(:,2),'.-','LineWidth',2)
% hold on
% plot(t(1) +3,t(2) + 1,'r*', 'MarkerSize',14);
% hold off
% grid on
% set(gca,'fontsize',16');
% xlabel('\nu_x');                 % Add labels
% ylabel('\nu_z');
% title(['Resonances up to 3rd order. Working Point ' num2str(t)])
% print('Tune_diagram','-dpng','-r300')

%%
% td1_x = [td1(1,1) td1(2,1)]
% td1_y = [td1(1,2) td1(2,2)]

figure

w = 1:2:length(td2)

for i=1:length(w)
qqq{i} = td2(w(i):w(i)+1,:)
r = qqq{i}
h = line(r(:,1),r(:,2),'LineWidth',2, 'Color','m');
hold on
end

w = 1:2:length(td1)

for i=1:4
qqq{i} = td1(w(i):w(i)+1,:)
r = qqq{i}
h=line(r(:,1),r(:,2),'LineWidth',2, 'Color','k');
h.LineStyle = '-'
hold on
end

plot(t(1) +3,t(2) + 1,'r*', 'MarkerSize',16);
plot(t(1) +3,t(2) + 1,'r.', 'MarkerSize',16);
plot(3.16,1.58,'b*', 'MarkerSize',16);
plot(3.16,1.58,'b.', 'MarkerSize',16);
hold off
grid on
set(gca,'fontsize',16');
xlabel('\nu_x');                 % Add labels
ylabel('\nu_z');
title(['Resonances up to 2rd order. Working Point ' num2str(t)])
print('Tune_diagram','-dpng','-r300')


%% Problem with 3rd order

figure
w = 1:2:length(td3)

figure
for i=1:2%length(w)-1
qqq{i} = td3(w(i):w(i)+1,:)
r = qqq{i}
line(r(:,1),r(:,2));
hold on
end

%%

% figure
% set(gcf,'color','w')
% line(td1_x,td1_y);
% hold on
% plot(td2(:,1), td2(:,2),'r.','LineWidth',2)
% %plot(td3(:,1), td3(:,2),'r.-','LineWidth',2)
% %plot(td4(:,1), td4(:,2),'g.-','LineWidth',2)
% %plot(td5(:,1), td5(:,2),'r.--','LineWidth',2)
% plot(t(1) +3,t(2) + 1,'r*', 'MarkerSize',14);
% hold off
% grid on
% set(gca,'fontsize',16');
% xlabel('\nu_x');                 % Add labels
% ylabel('\nu_z');
% title(['Resonances up to 3rd order. Working Point ' num2str(t)])
% %print('Tune_diagram','-dpng','-r300')

%% Not good => check!


R1 = tunespaceplot([3 4],[1 2],1);
R2 = tunespaceplot([3 4],[1 2],2);
R3 = tunespaceplot([3 4],[1 2],3);

figure

NL3 = size(R3,1);
for i = 1:NL3
    %hl = line([X1(i) X2(i)],[Y1(i) Y2(i)]);
    plot([R3(i,5) R3(i,7)],[R3(i,6) R3(i,8)],'m-','LineWidth',2);
  hold on
end

NL2 = size(R2,1);
for i = 1:NL2
    %hl = line([X1(i) X2(i)],[Y1(i) Y2(i)]);
    plot([R2(i,5) R2(i,7)],[R2(i,6) R2(i,8)],'b:','LineWidth',2);
  hold on
end

NL1 = size(R1,1);
for i = 1:NL1
    %hl = line([X1(i) X2(i)],[Y1(i) Y2(i)]);
    plot([R1(i,5) R1(i,7)],[R1(i,6) R1(i,8)],'k--','LineWidth',2.5);
   hold on
   
end
plot(t(1) +3,t(2) + 1,'r*', 'MarkerSize',16);
plot(t(1) +3,t(2) + 1,'r.', 'MarkerSize',16);
hold off
grid on
set(gca,'fontsize',16');
xlabel('\nu_x');                 % Add labels
ylabel('\nu_z');
title(['Resonances up to 3rd order. Working Point ' num2str(t)])
print('Tune_diagram','-dpng','-r300')
