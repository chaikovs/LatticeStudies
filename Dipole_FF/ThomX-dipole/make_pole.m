function [ pole ] = make_pole(radius,dx,teta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% make pole line
% input in m rad
% ouput in mm
teta=teta/2;
phi   =(-teta:teta/100:teta);
rr=radius+dx;
arc1   =[-rr*(1-cos(phi))+dx;rr*sin(phi)]; % X,S
rr=radius-dx;
arc2   =[-rr*(1-cos(-phi))-dx;rr*sin(-phi)]; % X,S
pole  =[arc1 arc2 arc1(:,1)]*1000;

end

