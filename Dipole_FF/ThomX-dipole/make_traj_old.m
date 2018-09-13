function [ traj, ds ] = make_traj(radius,dx,teta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% make half arc+straigth
% input in m rad
% ouput in mm

phi   =(0:teta/100:teta);
arc   =[-radius*(1-cos(phi))+dx;radius*sin(phi)]; % X,S
ds    =teta/100*radius;
str   =[arc(1,101)-ds*sin(teta)*(1:101) ; arc(2,101)+ds*cos(teta)*(1:101)];
traj  =[arc str]*1000;

end

