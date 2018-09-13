function get_H(RING)
%UNTITLED2 Summary of this function goes here
%   plot H function

[H]=CurlyH(RING,0,1:length(RING));
spos  = findspos(RING,1:length(RING));

figure(101);
set(gcf,'color','w')
plot(spos,H,'-g','LineWidth',2);
xlabel('S (m)')
ylabel('H (m)')
grid on

end

