function p = checkstability(fh_ring)

global THERING;

unstable_machines = 0;
p.nr_machines = 10;
THERING_error_free = THERING;
quaderr = 0:0.002:0.05;
p.quaderr = quaderr;%0:0.005:0.05;
dp = 0;

%% Loops over random machines

fprintf('Looping over random machines ...\n');
for ierr = 1:length(quaderr)
    fprintf('[The QUAD error is %1.3f]\n',quaderr(ierr));
    
    for i=1:p.nr_machines
        
        if rem(i,20) == 0
            fprintf('[Machine #%i/%i]\n', i, p.nr_machines);
        end
        
        
        % generate lattice errors
        %fprintf('Generating lattice errors ...\n');
        THERING = THERING_error_free;
        %getsp('QP1')
        seedquads(quaderr(ierr));
        % getsp('QP1')
        % tries calculation
        try
            
            %fprintf('Calculating matrix ...\n');
            
            [tx, ty] = tracem44(dp);
            if any(abs(tx)>2||abs(ty)>2)
                error('error in M44 matrix');
            end
        catch
            
            unstable_machines = unstable_machines + 1;
            %unstable_machine = true;
        end
        
        p.results.tx{ierr}{i} = tx;
        p.results.ty{ierr}{i} = ty;
        
    end
    
    p.unstable_machines{ierr} = unstable_machines;
    p.fraction_unstability{ierr} = unstable_machines/p.nr_machines;
    unstable_machines = 0;
    
end

%return;
save('stabilityReport.mat', 'p');
THERING = THERING_error_free;
fprintf('\n');
end
