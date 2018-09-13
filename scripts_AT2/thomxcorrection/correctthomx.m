% test errors and correction functions
close all
clear all

%addpath(genpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/correction/'));
%addpath(genpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/'));

% load lattice
ring=THOMX();

[l,t,c]=atlinopt(ring,0,1);

%% get RM
speclab='Chain';

% get indexes
indBPM=find(atgetcells(ring,'FamName','BPMx'))';
indHCor=find(atgetcells(ring,'FamName','HCOR'))';
indVCor=find(atgetcells(ring,'FamName','VCOR'))';
%indSCor=find(atgetcells(ring,'Class','Quadrupole'))';
indQCor=find(atgetcells(ring,'Class','Quadrupole'))';

ring=atsetfieldvalues(ring,indHCor,'PassMethod','StrMPoleSymplectic4Pass');
ring=atsetfieldvalues(ring,indVCor,'PassMethod','StrMPoleSymplectic4Pass');


modelrmfile=fullfile(pwd,['RMmodel' speclab '.mat']);%

if ~exist([modelrmfile],'file')
    
    ModelRM...
        =getresponsematrices(...
        ring,...
        indBPM,...
        indHCor,...
        indVCor,...
        indSCor,...
        indQCor,...
        [],...
        [0 0 0 0 0 0]',...
        [1 2 3 4 5 6 7 8 9 10 11 12]); % all RM
    
    save([modelrmfile],'ModelRM');
else
    load([modelrmfile],'ModelRM');
end


% mark quadrupoles to use for tune matching
indqf1=find(atgetcells(ring,'FamName','QF1\w*'));
ring=atsetfieldvalues(ring,indqf1,'ForTuneF',1);                
indqd2=find(atgetcells(ring,'FamName','QD2\w*'));
ring=atsetfieldvalues(ring,indqd2,'ForTuneD',1);                

inddq=find(atgetcells(ring,'FamName','DQ\w*'))';
inddl=find(atgetcells(ring,'FamName','DL\w*_3\w*'))';
ring=atsetfieldvalues(ring,[inddq inddl],'FitElement',1);     %mark as fitting point only some dipoles central ones.           

indDip=find(atgetcells(ring,'Class','Bend') & atgetcells(ring,'FitElement') )';

r0=ring;

zz=zeros(size(indBPM))';
r0=atsetbpmerr(r0,indBPM, zz,zz, zz,zz, zz,zz, zz);

% set errors, large, AT does not find a closed orbit
ind=find(atgetcells(ring,'Class','Quadrupole','Sextupole'));
dx=0.5e-3*randn(size(ind));
dy=0.5e-3*randn(size(ind));
dr=1.0e-3*randn(size(ind));

rerr=atsetshift(ring,ind,dx,dy);
rerr=atsettilt(rerr,ind,dr);
pb2q=atgetfieldvalues(rerr,indQCor,'PolynomB',{1,2});
rerr=atsetfieldvalues(rerr,indQCor,'PolynomB',{1,2},pb2q.*(1+randn(size(indQCor))*1e-2)');
%% correction chain

neigenvectors=[...
    10,... % n eig orbit H
    10,... % n eig orbit V
    15,... % skew quadrupole 
    15,... % normal quadrupole 
    10,... % fit normal quadrupole 
    10,... % fit dipole 
    10,... % fit skew quadrupole 
    ]; % number of eigenvectors 

diary('CorrChain.txt');
cororder=[1 2 3 6 1 2 3 6 -1];
%  '(-1 ): RF cavity frequency and time lag tuning '...
%  '( 0 ): open trajectory (finds closed orbit) '...
%  '( 1 ): orbit '...
%  '( 2 ): tune '...
%  '( 3 ): chromaticity '...
%  '( 4 ): dispersion '...
%  '( 5 ): dispersion free steering '...
%  '( 6 ): rdt + dispersion correction '...

rcor=CorrectionChain(...
    rerr,...            %1  initial lattice
    r0,...              %2  model lattice
    indBPM,...          %3  bpm index
    indHCor,...         %4  h steerers index
    indVCor,...         %5  v steerers index
    indSCor,...  %6  skew quad index
    indQCor,...      %7  quadrupole correctors index
    neigenvectors,...            %8  number of eigen vectors [NeigorbitH, NeigorbitV, NeigQuadrdt, Neigdispv, Neigdisph,neig rdt corr, SkewQuadRDT]
    cororder,...       %9  correction order 1: orbit, 2: tune, 3: skewquad disp v 4: quad disp h 5: quad RDT 6: skew RDT
    ModelRM,...          %10 response matrices
    '',...          %11 response matrices
    false);

diary off

%%
indBPM=find(atgetcells(ring,'FamName','BPMx'))';
indBPM=1:length(r0);

[l,t,c]=atlinopt(r0,0,indBPM);
[le,te,ce]=atlinopt(rerr,0,indBPM);
[lc,tc,cc]=atlinopt(rcor,0,indBPM);
bx0=arrayfun(@(a)a.beta(1),l);
bxe=arrayfun(@(a)a.beta(1),le);
bxc=arrayfun(@(a)a.beta(1),lc);
figure;
plot(bx0); hold on; 
plot(bxe); 
plot(bxc);
legend('initial','errors','corrected')
bbx=arrayfun(@(a,b)(a.beta(1)-b.beta(1))/b.beta(1),le,l);
bbxe=arrayfun(@(a,b)(a.beta(1)-b.beta(1))/b.beta(1),le,l);
bbxc=arrayfun(@(a,b)(a.beta(1)-b.beta(1))/b.beta(1),lc,l);
figure; plot(bbxe); hold on; plot(bbxc); legend('errors','corrected')


return