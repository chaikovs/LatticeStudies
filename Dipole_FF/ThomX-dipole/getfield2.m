function [B] = getfield2(X)
%GETFIELD Summary of this function goes here
%   Detailed explanation goes here
global  Xb Sb Bz 
% from OPERA field in 2D

% from file
bz=interp2(Xb,Sb,Bz,X(1)*1000,X(3)*1000);
%bs=interp1(S,Bsp,X(3))*X(2); % from gradient
B=[0 bz 0];