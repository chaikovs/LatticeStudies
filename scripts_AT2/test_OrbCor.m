

%%   Put in a vertical orbit distortion 

setsp('VCOR',0)    

Y0 = getam('BPMz'); 
BPMindex=getlist('BPMz');	

% Create an Orbit Error
vcm = .00005 * randn(12,1);  % ?????? Why 5 and not 12
setsp('VCOR', vcm);
% Get the vertical orbit
Y = getam('BPMz');	
figure
plot(Y-Y0);hold on;

%%

% Get the Vertical response matrix from the model
%Ry = measrespmat('BPMz', getlist('BPMz'), 'VCOR');  %getrespmat     % 13x12 matrix
Ry = getrespmat('BPMz', getlist('BPMz'), 'VCOR');  %getrespmat     % 13x12 matrix
% Computes the SVD of the response matrix
%%
Ivec = 1:10;
[U, S, V] = svd(Ry, 0);	


%%
% Find the corrector changes use the singular values
DeltaAmps = -V(:,Ivec) * S(Ivec,Ivec)^-1 * U(:,Ivec)' *  (Y-Y0);
% Changes the corrector strengths 
stepsp('VCOR', DeltaAmps);
%%
% Get the vertical orbit
Y = getam('BPMz');	
plot(Y-Y0, 'r');	