function [S, Bzn, Bsp] = getEXCELfield1D(file,brho,teta,Lbore)
%GETRADIAFIELD Summary of this function goes here
%   brho,teta to normalize field & deviation
%   remove the hard edge model


A=load(file);
Bz=A(:);
len=length(Bz);
S=(-685:5:685)/1000;
gap=0.038;

% % fit bore length
% n=22;% mm removed
% Bz(n+1:400)=Bz(1:400-n);Bz(1:n)=0;
% Bz=[Bz(1:400) ; Bz(401) ; Bz(400:-1:1) ]';

%
Bz=Bz/sum(Bz)*1000/5*(teta*brho);
% hard edge model
Bmax=teta/Lbore*brho % to get teta
rho=brho/Bmax

Bz0=0*Bz;
for i=1:len
    if (S(i)>=-Lbore/2) && (S(i)<=+Lbore/2)
        Bz0(i)=Bmax*300/300; % keep integral : sum # int
    end
end

%Bzn=Bz;
Bzn=Bz-Bz0; % remove hard edge model for tracking
Bsp=gradient(Bz)*1000/5; %Bs gradient
K= sum(Bz(138:275).*(Bz0(138)-Bz(138:275)))*0.005/Bz0(138)^2/gap
%
Kc=sum(Bz(138:239).*(Bz0(138)-Bz(138:239)))*0.005/Bz0(138)^2/gap
Kc=Kc-sum(Bz(239:275).*Bz(239:275))*0.005/Bz0(138)^2/gap
% figure(5)
% plot(S,Bsp*0.001)
% grid on


figure(1)
set(gca,'Fontsize',14)
plot(S*1000,Bz0,'-k');hold on
plot(S*1000,Bz,'-b');hold off
xlim([S(138) S(len)]*1000)
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


