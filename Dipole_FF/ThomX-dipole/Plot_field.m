function  Plot_field
%PLOT_FIELD Summary of this function goes here
%   Detailed explanation goes here
%   Maille ThomX

file1='Bprofil-s-ThomX-Dpole.dat';
fid = fopen(file1, 'r');
B_d = fscanf(fid, '%g %g', [2 inf]) ;
fclose(fid);

file1='Bprofil-s-ThomX-Qpole.dat';
fid = fopen(file1, 'r');
B_q = fscanf(fid, '%g %g', [2 inf]) ;
B_q(1,:)=B_q(1,:)-150;
fclose(fid);

file1='Bprofil-s-ThomX-Spole.dat';
fid = fopen(file1, 'r');
B_s = fscanf(fid, '%g %g', [2 inf]) ;
B_s(1,:)=B_s(1,:)-80;
fclose(fid);

figure(1)
stairs(B_q(1,:),B_q(2,:))
grid on


figure(2)
plot(B_d(1,:)-297.46/2,B_d(2,:),'r');hold on
plot(B_s(1,:)+140,B_s(2,:)*1000,'m');
plot(B_q(1,:)+275,B_q(2,:)*30,'b');
plot(B_s(1,:)+410,B_s(2,:)*1000,'m');
plot(B_q(1,:)+625,B_q(2,:)*30,'b');
plot(B_s(1,:)+550+150+115,B_s(2,:)*1000,'m');
plot(B_q(1,:)+1005,B_q(2,:)*30,'b');
plot(B_q(1,:)+1355,B_q(2,:)*30,'b');
plot(B_d(1,:)+1555+75+297.46/2,B_d(2,:),'r');

hold off
grid on

