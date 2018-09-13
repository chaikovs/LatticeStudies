function [S, Bzn, Bsp] = getMODELfield1D(file,brho,teta,Lbore)
%GETRADIAFIELD Summary of this function goes here
%   brho,teta to normalize field & deviation
%   remove the hard edge model

c=100;
S=(-400*c:1:400*c)/1000/c;
len=length(S);
gap=0.040;


% linear drop model
ng=1;
bb1=ones(1,(400-ng*gap*1000-(150-ng*gap*1000/2))*c)*0;
bb2=(1:ng*gap*1000*c)/(ng*gap*1000*c);
bb3=ones(1,(150-ng*gap*1000/2)*c);
bb=[bb1 bb2 bb3];
Bz=[bb 1 bb(400*c:-1:1)];
%length(Bz)

%
Bz=Bz/sum(Bz)*1000*c*(teta*brho);
%Bz(40000:70000)=Bz(42000:72000);

% hard edge model
Bmax=teta/Lbore*brho % to get teta
rho=brho/Bmax

Bz0=0*Bz;
for i=1:len
    if (S(i)>=-Lbore/2) && (S(i)<=+Lbore/2)
        Bz0(i)=Bmax*300/300.0; % keep integral : sum # int
    end
end


%Bzn=Bz;
Bzn=Bz-Bz0; % remove hard edge model for tracking
Bsp=gradient(Bz)*1000*c; %Bs gradient
K=sum(Bz.*(Bz(401*c)-Bz))*0.001/c/Bz(401*c)^2/gap/2

% figure(5)
% plot(S,Bsp*0.001)
% grid on


figure(1)
set(gca,'Fontsize',14)
plot(S*1000,Bz0,'-k');hold on
plot(S*1000,Bz,'-b');hold off
xlim([S(50000) S(60000)]*1000)
xlabel('s (mm)')
ylabel('B (T)')
legend('Hard edge profile','Real profile')
grid on

return

% K=sum(Bz.*(Bz(301)-Bz))*0.001/Bz(301)^2/0.04/2
% Lff=6*0.04*K
% 
% figure(3)
% plot(S,Bz.*(Bz(301)-Bz))
% grid on

% %
% ng=1;
% bb1=ones(1,400-ng*gap*1000-(150-ng*gap*500))*0;
% bb2=(1:ng*gap*1000)/(ng*gap*1000);
% bb3=ones(1,(150-ng*gap*500));
% bb=[bb1 bb2 bb3];
% Bz=[bb 1 bb(400:-1:1)];
% %length(Bz)


