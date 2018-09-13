function [s_alba,b2_alba,b3_alba,b4_alba,b5_alba] = getMEASUREDMultipolesSep(file)
% Read and plot Integrated field
% Fit ouput of  [Xb,Sb,Bz] = getOPERAfield2D(file)
sheet='Multipoles';  % 

H = xlsread(file,sheet);H=H'; % To match OPERA getfield
ss=size(H);
s_alba = H(1,:)*1e3; 
b2_alba = H(8,:);
b3_alba = H(10,:);
b4_alba = H(12,:);
b5_alba = H(14,:);
