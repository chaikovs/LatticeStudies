function [x,Bint,dBB] = getMEASUREDMultipoles(file)
% Read and plot Integrated field
% Fit ouput of  [Xb,Sb,Bz] = getOPERAfield2D(file)
sheet='DeltaB over B';  % 

H = xlsread(file,sheet);H=H'; % To match OPERA getfield
ss=size(H);
x = H(1,:)*1e3; 
dBB = H(4,:);
Bint = H(2,:);
