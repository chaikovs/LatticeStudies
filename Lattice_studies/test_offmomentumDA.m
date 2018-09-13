
%%

ring = ThomX_017_064_r56_02_chro00();

%%
% zzda= [];
% xxda= [];
% dp = -0.03:0.01:0.03;
% for idp = dp
% [xx,zz]=atdynap(RING,25,idp,0.02);
% 
% xxdai=[-10:0.2:10]*1e-3;
% zzdai=interp1(xx,zz,xxdai','pchip',0);
% 
% 
% zzda = [zzda zzdai(51)];
% 
% [indx, qq] = find(zz>=0.001);
% 
% if ~isempty(indx)
% xxda = [xxda xx(indx(1))];
% else
%     xxda = [xxda 0];
% end
% 
% end
% 
% %%
% 
% figure
% plot(1e2*dp,1e3*zzda,'.-' )
% ylim([0 30])
% 
% figure
% plot(1e2*dp,1e3*xxda,'r.-' )
% ylim([0 40])

%%
ring = RING;

dp = -0.05:0.01:0.05;
xlist = (0:2:38)*1e-3;
xmax = 0.0;
xoff = [];
for idp = dp
for i=1:length(xlist)
   rin=[xlist(i);0;1e-12;0;idp;0];
   [dummy,lost]=ringpass(ring,rin,500,'KeepLattice'); 
   if lost, break; end
end
if i>1 ; xmax=xlist(i-1);end
xoff = [xoff xmax];
end

%%
figure
plot(1e2*dp,1e3*xoff,'r.-' )
ylim([0 40])




