function Family = gethcmfamily
%GETHCMFAMILY - Returns the default horizontal corrector family
%  Family = gethcmfamily
%
%  The family name is determined by the MemberOf field equal to 'HCM' 
%
%  See also gethbpmfamily, getvbpmfamily, getvcmfamily
%
%  Writen by Greg Portmann

persistent WarningFlag 

Family = findmemberof('HCM');

if isempty(Family)
    Family = findmemberof('HCOR');
    if isempty(Family)
        %Family = {'HCM'};
        if isempty(WarningFlag)
            fprintf('\n   No default horizontal corrector family has been specified in the MML.\n');
            fprintf('   To define one, add ''HCM'' or ''HCOR'' to the .MemberOf field for the default family.\n\n');
            WarningFlag = 1;
        end

        Family = {''};
    end
end

Family = Family{1};

