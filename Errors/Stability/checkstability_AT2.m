function p = checkstability_AT2(fh_ring)

ring_woErr = ThomX_017_064_r56_02_chro00_AT2();
ring = ThomX_017_064_r56_02_chro00_multip_AT2();

unstable_machines = 0;
p.nr_machines = 20;
%THERING_error_free = THERING;
quaderr = 1e-12:0.002:0.05;
p.quaderr = quaderr;%0:0.005:0.05;
dp = 0;

r0 = ring;
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
        
        %getsp('QP1')
        rerr = seedquads_AT2(r0, quaderr(ierr));
        % getsp('QP1')
        % tries calculation
        try
            
            %fprintf('Calculating matrix ...\n');
            
            [tx, ty] = tracem44_AT2(rerr, dp);
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

fprintf('\n');
end
