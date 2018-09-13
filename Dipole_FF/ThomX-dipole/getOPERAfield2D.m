function [Xb,Sb,Bz] = getOPERAfield2D(file)
%GETRADIAFIELD Summary of this function goes here
%   Detailed explanation goes here

%file='field/champ_dipole_OPERA.table';

fid = fopen(file);
A = fscanf(fid, '%g %g %g %g', [4 inf]);
fclose(fid);

% n1=1001+1000;
% n2=201+200;

n1=1001;
n2=201;

Xb =reshape(A(1,:),n1,n2);
Sb =reshape(A(3,:),n1,n2);
Bz =reshape(A(4,:),n1,n2);
return

