% This script simulates preps through a simple phantom, then attempts to 
% reconstruct the phantom with the Cai2013 method, internally using a basis of
% synthetic materials to speed up convergence

% Set the size of the phantom
ObjectSize = 256;
MaxIters = 200;
NbPixelsPerProj = round(ObjectSize * sqrt(2));

% Simulate noisy preps
noise = true;
[y, k_d, A, M, S, T, ~] = SimulatePreps(ObjectSize, noise, true, 'fessler');
%%
S = repmat(S, [1 1 NbPixelsPerProj]); % Same spectrum, repeated for all pixels. Just to test for pixel-dependent spectrum feature

% Set regularization parameters
regulWeights = [1000000, 1000000, 100];
delta_h = [0.001, 0.001, 0.1];

% With regularization
doRegul = true;
runOnGPU = true;
useNesterov = false;
tic
[Cai2013_iterates, Cai2013_costs ]= Cai2013(y, ObjectSize, A, M, S, k_d, MaxIters, T, regulWeights, delta_h, runOnGPU, useNesterov);
toc

% Multiply the results by each material's density, in order to get something in g/ml
Cai2013_iterates(:,1,:)=Cai2013_iterates(:,1,:)*4.933; %Iodine
Cai2013_iterates(:,2,:)=Cai2013_iterates(:,2,:)*7.9; %Gadolinium
% Nothing to do for water, since it has density 1

% Show the resulting sequence of iterates
PlayIterates(Cai2013_iterates, [0 0 0], [0.015 0.015 1.5])

% Plot the resulting costs, in loglog scale
loglog(Cai2013_costs - min(Cai2013_costs(:)));

% % Save the iterates, if necessary
% save('../Results/Cai2013NoRegul.mat', 'Cai2013_noRegul_iterates');
% save('../Results/Cai2013Regul.mat', 'Cai2013_Regul_iterates');
