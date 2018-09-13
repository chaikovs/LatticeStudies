close all;
clear all;
fid = fopen('Classeur1.txt');
%for i=1:3
%ligne=fgetl(fid);
%end
a = fscanf(fid,'%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t ',[11 Inf]);
a=a';
MeanDCCT=a(:,1);
%num : numero dipole
B=a(:,2);
%B0 : champ au centre 
IntB1=a(:,3);
%Angle : Angle de déviation
IntA1=a(:,4);
IntB2=a(:,5);
IntA2=a(:,6);
IntB3=a(:,7);
IntA3=a(:,8);
IntB4=a(:,9);
IntA4=a(:,10);
Roll=a(:,11);

Il=a(1:14,1);
Byl=a(1:14,2);

figure
    plot(MeanDCCT,Roll,'+');
    axis([0,300,-0.3,0.2]) %Axes ajustés
    grid on
    title('Roll Angle en fonction du courant');
 	xlabel ('Courant (A)');ylabel('Angle de défaut (°)');
    set(gca,'fontsize',14)
    %set(title,'fontsize',18)
    legend('Roll Angle Dipole#09')
    hold on
    
 figure
    plot(MeanDCCT,IntA1,'+');
    axis([0,300,-9e-4,2e-4]) %Axes ajustés
    grid on
    title('Intégrale de A1 en fonction du courant');
 	xlabel ('Courant (A)');ylabel('Int. A1 (T.m)');
    set(gca,'fontsize',14)
    %set(title,'fontsize',18)
    legend('Int A1 Dipole#09')
    hold on
    
     figure
    plot(MeanDCCT,IntA2,'+');
    axis([0,310,-1e-4,6e-4]) %Axes ajustés
    grid on
    title('Intégrale de A2 en fonction du courant');
 	xlabel ('Courant (A)');ylabel('Int. A2 (T.m)');
    set(gca,'fontsize',14)
    %set(title,'fontsize',18)
    legend('Int A2 Dipole#09')
    hold on
    
   figure
    plot(MeanDCCT,IntA3,'+');
    axis([0,310,-6e-5,2e-5]) %Axes ajustés
    grid on
    title('Intégrale de A3 en fonction du courant');
 	xlabel ('Courant (A)');ylabel('Int. A3 (T.m)');
    set(gca,'fontsize',14)
    %set(title,'fontsize',18)
    legend('Int A3 Dipole#09')
    hold on
    
    
   figure
    plot(MeanDCCT,IntA4,'+');
    axis([0,310,-1e-5,8e-5]) %Axes ajustés
    grid on
    title('Intégrale de A4 en fonction du courant');
 	xlabel ('Courant (A)');ylabel('Int. A4 (T.m)');
    set(gca,'fontsize',14)
    %set(title,'fontsize',18)
    legend('Int A4 Dipole#09')
    hold on 
    
     figure
    plot(MeanDCCT,IntB1,'+');
    axis([0,300,0,0.25]) %Axes ajustés
    grid on
    title('Intégrale de B1 en fonction du courant');
 	xlabel ('Courant (A)');ylabel('Int. B1 (T.m)');
    set(gca,'fontsize',14)
    %set(title,'fontsize',18)
    legend('Int B1 Dipole#09')
    hold on
    
    figure
    plot(MeanDCCT,IntB2,'+');
    axis([0,300,-4e-4,1e-4]) %Axes ajustés
    grid on
    title('Intégrale de B2 en fonction du courant');
 	xlabel ('Courant (A)');ylabel('Int. B2 (T.m)');
    set(gca,'fontsize',14)
    %set(title,'fontsize',18)
    legend('Int B2 Dipole#09')
    hold on
    
    figure
    plot(MeanDCCT,IntB3,'+');
    axis([0,300,-4e-4,1e-4]) %Axes ajustés
    grid on
    title('Intégrale de B3 en fonction du courant');
 	xlabel ('Courant (A)');ylabel('Int. B3 (T.m)');
    set(gca,'fontsize',14)
    %set(title,'fontsize',18)
    legend('Int B3 Dipole#09')
    hold on
    
    figure
    plot(MeanDCCT,IntB4,'+');
    axis([0,300,-3e-5,0]) %Axes ajustés
    grid on
    title('Intégrale de B4 en fonction du courant');
 	xlabel ('Courant (A)');ylabel('Int. B4 (T.m)');
    set(gca,'fontsize',14)
    %set(title,'fontsize',18)
    legend('Int B4 Dipole#09')
    hold on
     
%     plot(MeanDCCT,B,'b');
%     hold on
%     plot(Il,Byl,'r');
%     %plot (Il,Byl,'r');
%     axis tight %Axes ajustés
%     grid on
%     title('Courbe excitation By(I)');
%     set(gca,'fontsize',14)
%  	xlabel ('I(A)');ylabel('By(G)');
%     legend('By de 0A à 300A','By lineaire de 0A à 140A','Fit Lmag +/- 10mm','fontsize',14)
%     hold off