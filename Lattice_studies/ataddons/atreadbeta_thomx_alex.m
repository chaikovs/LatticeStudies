function [superp,periods]=atreadbeta_thomx_alex(filename,cavipass,bendpass,quadpass)
%ATREADBETA				reads a BETA file
%
%ring=ATREADBETA(fname,cavipass,bendpass,quadpass,multipass)
%
%FILENAME:	BETA file
%CAVIPASS:	pass method for cavities (default IdentityPass)
%BENDPASS:	pass method for dipoles (default BndMPoleSymplectic4E2Pass)
%QUADPASS:	pass method for quadrupoles (default QuadMPoleFringePass)
%MULTIPASS:	pass method for sextupoles (default StrMPoleSymplectic4Pass)
%
%[superp,periods]=ATREADBETA(fname,cavipass,bendpass,quadpass)
%		returns only one superperiod and the number of superperiods
%
%See also: ATX, ATLINOPT, ATMODUL

global GLOBVAL
persistent fpath

if isempty(fpath), fpath=getenv('DBETA'); end

if nargin < 5, multipass='StrMPoleSymplectic4Pass'; end
if nargin < 4, quadpass= 'StrMPoleSymplectic4Pass'; end
if nargin < 3, bendpass= 'BndMPoleSymplectic4Pass'; end
%if nargin < 3, bendpass= 'BndMPoleSymplecticNew4Pass'; end
if nargin < 2, cavipass= 'IdentityPass'; end
if nargin < 1, filename=''; end

if isempty(filename)
    [fname,fpath]=uigetfile('*.str','BETA structure',[fpath filesep]);
    if ~ischar(fname), error('ReadBeta:NoFile','No file selected'); end
    filename=fullfile(fpath,fname);
end

fid=fopen(filename,'rt');
if fid < 0, error('ReadBeta:NoFile','Cannot open file %s',filename); end

betadelim(fid,'LIST OF ELEMENTS');
GLOBVAL.E0=50E6;
line=fgetl(fid);
nb_elems=sscanf(line,'%d');
elemlist=cell(nb_elems,1);
eldict=cell(nb_elems,1);
for el=1:nb_elems
    nextelem=readelem(fid,cavipass,bendpass,quadpass,multipass);
   elemlist{el}=nextelem;
   eldict{el}=nextelem.FamName;
end
disp(['Elements processed (' num2str(nb_elems) ' elements)']);

betadelim(fid,'STRUCTURE');
line=fgetl(fid);
nb_stru=sscanf(line,'%d');
id_stru=0;
superp=cell(nb_stru+1,1);	% add 1 for implicit cavity
dipelem=[];
anglein=0;
angleout=0;
lff=0;
displ=zeros(1,3);
srot=0;
for el=1:nb_stru
   elcode=fscanf(fid,'%s',1);
   elnum=str2double(elcode);
   if ~isfinite(elnum)
      elnum=find(strcmp(elcode,eldict));   end
   try
      nextelem=elemlist{elnum};
   catch
      error('ReadBeta:BadElem',['Cannot identify element ' elcode]);
   end
   switch nextelem.BetaCode
   case 'CO'
      if isempty(dipelem)
	 anglein=nextelem.Angle;
	 lff=nextelem.Lff;
      else
	 angleout=nextelem.Angle;	% create immediately in case of 2 adjacent CO elems
         id_stru=id_stru+1;
	 superp{id_stru}=tunedipole(dipelem,anglein,angleout,lff);
	 anglein=0;
	 angleout=0;
	 lff=0;
	 dipelem=[];
      end
   case 'RO'
      srot=srot+nextelem.Srot;
   case 'DE'
      displ=displ+nextelem.Displacement;
   otherwise
      if ~isempty(dipelem)
         id_stru=id_stru+1;
	 superp{id_stru}=tunedipole(dipelem,anglein,angleout,lff);
	 anglein=0;
	 angleout=0;
	 lff=0;
	 dipelem=[];
      end
      if srot ~= 0
      if isfield(nextelem,'R1')
	 srollmat=mkSRotationMatrix(srot);
	 nextelem.R1=srollmat;
	 nextelem.R2=srollmat';
      end
      end
      if max(abs(displ)) > 0
      if isfield(nextelem,'T1')
         nextelem.T1([1 3 5])=displ;
         nextelem.T2([1 3 5])=-displ;
      end
      end
      if strcmp(nextelem.BetaCode,'DI')
	 dipelem=nextelem;
      else
	 id_stru=id_stru+1;
	 superp{id_stru}=nextelem;
      end
   end
end
superp(id_stru+1:end)=[];
disp(['structure processed (' num2str(length(superp)) ' elements)']);

nper=fscanf(fid,'%d',1);
fclose(fid);

% cavities=findcells(superp,'BetaCode','CA');
% 
% if isempty(cavities)		% add implicit cavity if necessary
%    famcav=[];
%    for i=1:length(elemlist)
%       if strcmp(elemlist{i}.BetaCode,'CA')
%         famcav=[famcav;i];
%       end
%    end			% add implicit RF cavity
%    if isempty(famcav)
%       famcav=rfcavity('RFCAV',0,0,0,1,cavipass);
%       elemlist{famcav}.BetaCode='CA';
%    end
%    superp{end+1}=elemlist{famcav(1)};
%    cavities=length(superp);
% end

superp=setcellstruct(superp,'Energy',1:length(superp),GLOBVAL.E0);

% if nargout >= 2
%    periods=nper;
%     for i=cavities			% set cavity frequency
%        superp{i}=tunecavity(superp{i},9E6/16/length(cavities),...
%         findspos(superp,id_stru+1),1,nper);
%     end
% else
%     for i=cavities			% set cavity frequency
%        superp{i}=tunecavity(superp{i},9E6/16/length(cavities),...
%         findspos(superp,id_stru+1),nper,nper);
%     end
%     if nper > 1
%        superp=repmat(superp,1,nper);
%     end
% end


evalin('base','global GLOBVAL');

function cav=tunecavity(cav0,V,clength,ncell,ntot)
cav=cav0;
frev=2.997924E8/clength;
if cav.Voltage == 0, cav.Voltage=V; end
if cav.HarmNumber > 1
    harm=ceil(cav.HarmNumber/ntot);
else
    harm=round(1.66678*clength);
end
cav.Frequency=frev*harm;
cav.HarmNumber=ncell*harm;

function newelem=readelem(fid,cavipass,bendpass,quadpass,multipass)

global GLOBVAL

line=fgetl(fid);
next=1;
[elname,count,errmess,nl]=sscanf(line(next:end),'%s',1);
next=next+nl;
[code,count,errmess,nl]=sscanf(line(next:end),'%s',1);
next=next+nl;
params=sscanf(line(next:end),'%f')';
switch (code)
case 'SD'
   newelem=atelem('drift','FamName',elname,'Length',params(1));
case 'QP'
   newelem=atelem('quadrupole','FamName',elname,'Length',params(1),...
   'K',params(2),'PolynomA',[0 0 0],'PolynomB',[0 params(2) 0],'MaxOrder',2,'NumIntSteps',10);
   newelem.PassMethod=quadpass;
case 'DE'
   newelem=atelem('marker','FamName',elname,'Displacement',params(1:3));
case 'RO'
   newelem=atelem('marker','FamName',elname,'Srot',params(1));
case 'CO'
   newelem=atelem('marker','FamName',elname,'Angle',params(1),'Lff',params(3));
case 'DI'
   strength=-params(3)/params(2)/params(2);
   newelem=atelem('rbend3','FamName',elname,'Length',params(1)*params(2),'BendingAngle',params(1),...
   'K',strength,'PolynomA',[0 0 0],'PolynomB',[0 strength 0],'FringeInt1',0.5,'FringeInt2',0.5, ...
   'MaxOrder',3,'NumIntSteps',10);
   newelem.PassMethod=bendpass;
case 'SX'
   if params(1) < 1e-10
      code='LD3';
      newelem=atelem('sextupole','FamName',elname,'PolynomB',[0 0 params(1)*params(2)],'PassMethod','ThinMPolePass');
   else
      newelem=atelem('sextupole','FamName',elname,'Length',params(1),'PolynomB',[0 0 params(2)],'PassMethod',multipass);
   end
case 'LD'
   order=params(2)/2;
   polb=[];
   polb(1,order)=params(1);
   code=[code int2str(order)];
   newelem=atelem('marker','FamName',elname,'PolynomA',zeros(1,order),'PolynomB',polb,'MaxOrder',order-1);
   newelem.PassMethod='ThinMPolePass';
case 'LT'
   order=params(2)/2;
   pola=[];
   pola(1,order)=params(1);
   code=[code int2str(order)];
   newelem=atelem('marker','FamName',elname,'PolynomA',pola,'PolynomB',zeros(1,order),'MaxOrder',order-1);
   newelem.PassMethod='ThinMPolePass';
case 'KI'
   if params(3) > 0, code='CHV'; end
   newelem=atelem('corrector','FamName',elname,'KickAngle',[params(1) params(2)]);
   newelem.PassMethod='IdentityPass';
case 'CA'
   GLOBVAL.E0=params(3);
   newelem=atelem('marker','FamName',elname,'Voltage',abs(params(1)),'Frequency',0,'HarmNumber',params(2));
   newelem.PassMethod=cavipass;
otherwise
   newelem=atelem('marker','FamName',elname);
end
newelem.BetaCode=code;

function betadelim(fid,code)
while true
   while true
      line=fgetl(fid);
      if ~ischar(line), error('ReadBeta:EndOfFile','Encountered unexpected end of file.'); end
      if ~isempty(strfind(line,'***')), break; end
   end
   if ~isempty(strfind(line,code)), break; end
end

function dip2=tunedipole(dip1,anglein,angleout,lff)
dip2=dip1;
dip2.EntranceAngle=anglein;
dip2.ExitAngle=angleout;
%dip2.FringeInt=lff;
dip2.FullGap=lff/3;
dip2.EdgeEffect1=1;
dip2.EdgeEffect2=1;
dip2.FringeBendEntrance = 3;
dip2.FringeBendExit = 3;


