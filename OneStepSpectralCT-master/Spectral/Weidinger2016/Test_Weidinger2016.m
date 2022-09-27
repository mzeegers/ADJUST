% This script simulates preps through a simple phantom, then attempts to 
% reconstruct the phantom with the Weidinger2016 method, internally using a basis of
% synthetic materials to speed up convergence

% Set the size of the phantom
ObjectSize = 64;
NbIters = 20;
NbPixelsPerProj = round(ObjectSize * sqrt(2));

% Simulate noisy preps
noise = true;
[y, ~, A, M, S, T, ~] = SimulatePreps_more_bins(ObjectSize, noise, false, 'none');
S = repmat(S, [1 1 NbPixelsPerProj]); % Same spectrum, repeated for all pixels. Just to test for pixel-dependent spectrum feature

% Set regularization parameters
lambda = [1, 1, 1];

% Perform reconstruction
runOnGPU = false;
size(y), ObjectSize, size(A), size(M), size(S), size(T), NbIters, lambda

[Weidinger2016_iterates, costs ]= Weidinger2016(y, ObjectSize, A, M, S, T, NbIters, lambda, runOnGPU);

% Multiply the results by each material's density, in order to get something in g/ml
Weidinger2016_iterates(:,1,:)=Weidinger2016_iterates(:,1,:)*4.933; %Iodine
Weidinger2016_iterates(:,2,:)=Weidinger2016_iterates(:,2,:)*7.9; %Gadolinium
% Nothing to do for water, since it has density 1

% Show the resulting sequence of iterates
PlayIterates(Weidinger2016_iterates, [0 0 0], [0.015 0.015 1.5])

% Plot the resulting costs, in loglog scale
loglog(costs - min(costs(:)))
