function [B] = getfield(X)
%GETFIELD Summary of this function goes here
%   Detailed explanation goes here
global  S Bz Bsp
% simple dipole
%B=[0 0.4411 0];


% from file
bz=interp1(S,Bz,X(3));
bs=interp1(S,Bsp,X(3))*X(2); % from gradient
B=[0 bz bs];