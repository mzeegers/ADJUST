% This script simulates preps through a simple phantom, then attempts to 
% reconstruct the phantom with the Barber2016 method

% Set the size of the phantom
ObjectSize = 32;
NbIters = 200;
NbPixelsPerProj = round(ObjectSize * sqrt(2));

% Simulate noisy preps
noise = true;
normalizeSpectrum = false;
[y, ~, A, M, S, T, ~, x] = SimulatePreps(ObjectSize, noise, normalizeSpectrum, 'normalize');
S = repmat(S, [1 1 NbPixelsPerProj]); % Same spectrum, repeated for all pixels. Just to test for pixel-dependent spectrum feature

% Regularization parameters (one per material)
gammas = [64 64 6400];

% Run Barber2016 with magical parameter lambda
lambda = 0.00003;
theta = 0.6;
tic
runOnGPU = true;
[Barber2016_iterates, gaps, costs]= Barber2016(y, ObjectSize, A, M, S, T, NbIters, lambda, theta, gammas, runOnGPU);
toc

% Multiply the results by each material's density, in order to get something in g/ml
Barber2016_iterates(:,1,:)=Barber2016_iterates(:,1,:)*4.933; %Iodine
Barber2016_iterates(:,2,:)=Barber2016_iterates(:,2,:)*7.9; %Gadolinium
% Nothing to do for water, since it has density 1

% Show the resulting sequence of iterates
PlayIterates(Barber2016_iterates, [0 0 0], [0.015 0.015 1.5])

% Plot the convergence curve
figure(1); loglog(abs(gaps));

% Plot the convergence curve
figure(2); loglog(costs - min(costs));

% % Save the iterates, if necessary
save('/home/mory/data/MatlabOneStep/Barber/Barber_64_6400.mat', 'Barber2016_iterates', '-v7.3');
