close all;
clear all;
fid = fopen('Classeur3.txt');
%for i=1:4
%ligne=fgetl(fid);
%end
a = fscanf(fid,'%f\t %f\t %f\t %f\t',[4 Inf]);
a=a';
num=a(:,1);
%num : numero dipole
B0=a(:,2);
%B0 : champ au centre 200A
IntB=a(:,3);
%IntB : Integrale de champ 200A
Lmag=a(:,4);
%Lmag : Longueur magnetique 200A  
%I=a(6:16,1);
%Lmag1=a(6:16,2);
figure
    plot(num,Lmag,'+');
    axis([0,16,295,296]) %Axes ajustés
    grid on
    title('Longueur magnétique à 200A');
 	xlabel ('N° Dipole');ylabel('Lmag(mm)');
    set(gca,'fontsize',14)
    %set(title,'fontsize',18)
    legend('Lmag')
    hold on
    
 figure
    plot(num,IntB,'+');
    axis([0,16,0.1574,0.1586]) %Axes ajustés
    grid on
    title('Intégrale de champ à 200A');
 	xlabel ('N° Dipole');ylabel('Int. B(T.m)');
    set(gca,'fontsize',14)
    %set(title,'fontsize',18)
    legend('Int. B')