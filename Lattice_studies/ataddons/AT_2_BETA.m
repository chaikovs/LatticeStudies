function AT_2_BETA(AT_ring,linename)
% function AT_2_OPA(AT_ring,linename)
% this functions converts the AT lattice AT_ring in OPA format.
% 
% 
% file ['' linename '_lattice.opa'] is generated contiaining the lattice
% elements definitions and the LINE. no other comands introduced
% 
%
% OPA may be found here: http://people.web.psi.ch/streun/opa/
%

outfile=['' linename '_lattice.str'];

%outfile='madXelemdef.elem';
%elelat=['{com   madX lattice elements: ' linename ' com}\n{com   Created: ' datestr(now) ' com}\n'];

%% get family names for definitions
[families,ind_first_oc_ring]=...
    unique(getcellstruct(AT_ring,'FamName',1:length(AT_ring)),'first');


elelat=['*** LIST OF ELEMENTS *** \n'  num2str(length(families)) '\n'  ];

format='%8.10f';

%% loop families for definitions
for i=1:length(families)
   el= AT_ring{ind_first_oc_ring(i)};
   if isfield(el,'BetaCode')
       type=el.BetaCode;
   elseif isfield(el,'Class')
       type=el.Class;
   else
       type='Marker';
   end
      
    switch type
        
        case {'DI','Dipole','Bend'} % dipole
            
        R=el.('Length')/el.('BendingAngle');
        ind=-el.('PolynomB')(2)*R*R;
        
            di=[' ' el.('FamName')   ' '...
                ' DI  ' num2str(el.('BendingAngle'),format)  ' '...
                num2str(R,format)  ' '...
                num2str(ind,format)  ' '...
               ];
           
            elelat=[elelat di '\n']; %#ok<*AGROW>
        case {'QP','Quadrupole'} % quadrupole
            
 
            qp=[' ' el.('FamName') '  '...
                ' QP  ' num2str(el.('Length'),format)  ' '...
                 num2str(el.('PolynomB')(2),format) ...
               ];
 
           
            elelat=[elelat qp '\n'];
        case {'SX','Sextupole'} % sextupole
            
            sx=[' ' el.('FamName') '  '...
                ' SX  ' num2str(el.('Length'),format)  ' '...
                 num2str(el.('PolynomB')(3),format) ...
               ];
           
            elelat=[elelat sx '\n'];
            
        case {'OC','Octupole'} % sextupole

        case {'MP','Multipole'} % multipole
           
        case {'ThinMultipole'} % multipole
        warning('still to be defined!')
        

        case {'DR','Drift','SD'} % drift
            
            dr=[' ' el.('FamName') '  '...
                ' SD  ' num2str(el.('Length'),format)  ' '...
               ];
           
            elelat=[elelat dr '\n']; 	
            
           
        otherwise
            warning(['Element: ' el.('FamName') ' was not converted, since it does not match any Class.'])
            mrk=[' ' el.('FamName') ': marker' '; '...
                ];
            elelat=[elelat mrk '\n\n'];
    end
  
end


%elelat=[elelat '\n\n {com ----- table of segments -------------- com}\n\n'];

elelat=[elelat '*** STRUCTURE *** \n'  num2str(length(AT_ring)) '\n'  ];

%% define lattice line
% loop all elements

n6=floor(length(AT_ring)/6)*6;
for i=1:6:n6

        elelat=[elelat '' AT_ring{i}.('FamName')  '  ' AT_ring{i+1}.('FamName')  '  '   ... 
            AT_ring{i+2}.('FamName')  '  ' AT_ring{i+3}.('FamName') '  '  ... 
            AT_ring{i+4}.('FamName')  '  ' AT_ring{i+5}.('FamName')   ... 
            '\n'];
        
end
for i=n6+1:length(AT_ring)

        elelat=[elelat '' AT_ring{i}.('FamName')  '\n'];
        
end



elelat=strrep(elelat,'RFC','RFCav');
%% print to file

of=fopen(outfile,'w');
fprintf(of,elelat);

fclose('all');




return