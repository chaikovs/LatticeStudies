close all;
clear all;
fid = fopen('Dipole09.txt');
for i=1:2
ligne=fgetl(fid);
end
a = fscanf(fid,'%f\t  %f\n',[2 Inf]);
a=a';
I=a(:,1);
By=a(:,2);
b = fscanf(fid,'%f\t  %f\n',[2 10]);
b=b';
Il=a(1:10,1);
Byl=a(1:10,2);
%Pm=polyfit(I,By,m)
fit1 = polyfit(a(1:10,1), a(1:10,2), 1)
figure

    %plot(I,By,'b');
    %hold on
    %plot(Il,Byl,'r');
    %plot (Il,Byl,'r');
    
    
    figure(1) 
    plot(I,By,'b');
    figure(2)  
    plot(Il,Byl,'r');
    
    axis tight %Axes ajustés
    grid on
    title('Courbe excitation By(I)');
 	xlabel ('I(A)');ylabel('By(G)');
    hold off