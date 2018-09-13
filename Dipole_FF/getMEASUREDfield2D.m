function [Xb,Sb,Bz] = getMEASUREDfield2D(file)
% Read and plot Integrated field
% Fit ouput of  [Xb,Sb,Bz] = getOPERAfield2D(file)
sheet='By';  % 

H = xlsread(file,sheet);H=H'; % To match OPERA getfield
ss=size(H);
Sb = H(:,1)*1e3; Sb(1)=[];   % over 301 point step= 2 mm
Xb = H(1,:)*1e3; Xb(1)=[];
[Xb, Sb]=meshgrid(Xb,Sb);
Bz = H((2:ss(1)),(2:ss(2)));


% Keep half of the map along S
% Sb=Sb((151:ss(2)-1));
% Bz=Bz(:,(151:ss(2)-1))';
