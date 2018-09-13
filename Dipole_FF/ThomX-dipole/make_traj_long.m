function [ traj, ds , np] = make_traj_long(radius,dx,teta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% make half arc+straigth
% input in m rad
% ouput in mm

np=200;
np1=np+1;


%
phi   =(0:teta/np:teta);
arc   =[-radius*(1-cos(phi))+dx;radius*sin(phi)]; % X,S
ds    =teta/np*radius*ones(1,np1);

%
ds0   = teta/np*(radius-dx)*(1:220);
%ds0   = teta/np*(radius)*(1:np1);
str   =[arc(1,np1)-ds0*sin(teta) ; arc(2,np1)+ds0*cos(teta)];
traj  =[arc str]*1000;

%
ds=[ds teta/np*(radius-dx)*ones(1,220)];
%ds=[ds teta/np*(radius)*ones(1,np1)];
end

