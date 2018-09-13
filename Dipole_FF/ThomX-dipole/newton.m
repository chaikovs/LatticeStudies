function [dX] = newton(t,X)
%NEWTON Summary of this function goes here
%   F=ma in magnetic field
% s step vector
% v input speed
% dv dv/ds
global E
A=1.7563e+11;  % q/m
gam=(E+0.551)/0.511;
A=A/gam;

% X coordinates
dX(1)=X(4);
dX(2)=X(5);
dX(3)=X(6);

[B] = getfield2(X(1:3));

% Speed
dX(4)=A*(X(5)*B(3)-X(6)*B(2));
dX(5)=A*(X(6)*B(1)-X(4)*B(3));
dX(6)=A*(X(4)*B(2)-X(5)*B(1));

dX=dX';

return