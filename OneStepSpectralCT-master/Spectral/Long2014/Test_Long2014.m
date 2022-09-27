% This script simulates preps through a simple phantom, then attempts to 
% reconstruct the phantom with the Long2014 method

% Set the size of the phantom
ObjectSize = 64;
NbIters = 10;
NbSubsets = 4;
NbPixelsPerProj = round(ObjectSize * sqrt(2));

% Simulate noisy preps
noise = true;
[y, ~, A, M, S, T, ~] = SimulatePreps_more_bins(ObjectSize, noise, false, 'none');
S = repmat(S, [1 1 NbPixelsPerProj]); % Same spectrum, repeated for all pixels. Just to test for pixel-dependent spectrum feature

% Set regularization parameters
lambda = [100, 100, 100];
delta_hyperbola = [0.0001, 0.0001, 0.01];

% - the LongAndFessler2014 method
runOnGPU = false;
tic
[Long2014_iterates, costs] = Long2014(y, ObjectSize, NbPixelsPerProj, A, M, S, T, NbIters, lambda, NbSubsets, delta_hyperbola, runOnGPU);
toc

% Multiply the results by each material's density, in order to get something in g/ml
Long2014_iterates(:,1,:)=Long2014_iterates(:,1,:)*4.933; %Iodine
Long2014_iterates(:,2,:)=Long2014_iterates(:,2,:)*7.9; %Gadolinium
% Nothing to do for water, since it has density 1

% Show the resulting sequence of iterates
PlayIterates(Long2014_iterates, [0 0 0], [0.015 0.015 1.5])

% Plot the resulting costs, in loglog scale
loglog(costs - min(costs(:)))

