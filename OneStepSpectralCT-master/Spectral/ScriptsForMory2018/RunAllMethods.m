% This script simulates preps through a simple phantom, then attempts to 
% reconstruct the phantom with five different one-step inversion methods:
% - Cai2013
% - Long2014
% - Weidinger2016
% - Mechlem2017
% - Barber2016

% Set the size of the phantom
ObjectSize = 256;
NbPixelsPerProj = round(ObjectSize * sqrt(2));

% Set regularization parameters
delta_huber = [0.001, 0.001, 0.1];
runOnGPU = true;

% Set the number of iterations for most methods
NbIters = 5000;

%% Cai2013
useNesterov = false;
% Load preps and repmat the spectrum
for mu_precond = {'none','normalize', 'orthonormalize'}
    
    load(['/home/mory/data/MatlabOneStep/Cai/', mu_precond{1},'_preps.mat'],'y','k_d', 'A', 'M', 'S', 'T');
    S = repmat(S, [1 1 NbPixelsPerProj]);
    
    % Run Cai2013
    regulWeights = [100000, 100000, 100];
    [Cai2013_iterates, Cai2013_costs ]= Cai2013(y, ObjectSize, A, M, S, k_d, NbIters, T, regulWeights, delta_huber, runOnGPU, useNesterov);

    % Save the results
    filename = ['/home/mory/data/MatlabOneStep/Cai/Cai2013_', mu_precond{1}, '_', num2str(regulWeights(1)), '_', num2str(regulWeights(2)), '_' , num2str(regulWeights(3)), '.mat'];
    save(filename, 'Cai2013_iterates', 'Cai2013_costs', '-v7.3');
end


%% SQS methods

% Load preps and repmat the spectrum
load('/home/mory/data/MatlabOneStep/Long/preps.mat', 'y', 'A', 'M', 'S', 'T');
S = repmat(S, [1 1 NbPixelsPerProj]);

% Run Weidinger2016
regulWeights = [10000, 10000, 10];
[Weidinger2016_iterates, Weidinger2016_costs ]= Weidinger2016(y, ObjectSize, A, M, S, T, NbIters, regulWeights, runOnGPU);

% Save the results
filename = ['/home/mory/data/MatlabOneStep/Weidinger/Weidinger2016_', num2str(regulWeights(1)), '_', num2str(regulWeights(2)), '_' , num2str(regulWeights(3)), '.mat'];
save(filename, 'Weidinger2016_iterates', 'Weidinger2016_costs', '-v7.3');

% Run Long2014
NbSubsets = 20; % For ordered subsets
regulWeights = [100000, 100000, 100];
[Long2014_iterates, Long2014_costs ]= Long2014(y, ObjectSize, NbPixelsPerProj, A, M, S, T, NbIters, regulWeights, NbSubsets, delta_huber, runOnGPU);

% Save the results
filename = ['/home/mory/data/MatlabOneStep/Long/Long2014_', num2str(regulWeights(1)), '_', num2str(regulWeights(2)), '_' , num2str(regulWeights(3)), '.mat'];
save(filename, 'Long2014_iterates', 'Long2014_costs', '-v7.3');

% Run Mechlem2017 with less iterations
NbIters = 200;
NbSubsets = 4; % More than that speeds up computation, but makes it unstable, and ultimately, divergent
regulWeights = [100000, 100000, 100]; % Huber is scaled by delta, so where there is a low delta, there should be a high lambda to compensate
[Mechlem2017_iterates, Mechlem2017_costs ]= Mechlem2017(y, ObjectSize, NbPixelsPerProj, A, M, S, T, NbIters, regulWeights, NbSubsets, delta_huber, runOnGPU);

% Save the results
filename = ['/home/mory/data/MatlabOneStep/Mechlem/Mechlem2017_', num2str(regulWeights(1)), '_', num2str(regulWeights(2)), '_' , num2str(regulWeights(3)), '.mat'];
save(filename, 'Mechlem2017_iterates', 'Mechlem2017_costs', '-v7.3');

%% Barber2016

% Simulate noisy preps
load('/home/mory/data/MatlabOneStep/Barber/preps.mat', 'y', 'A', 'M', 'S', 'T');

% Run Barber2016
NbIters = 10000;
lambda = 0.0001; % Controls convergence
theta = 0.5; % Between 0 and 1, see Chambolle Pock
TVthresholds = [100,100,10000];
[Barber2016_iterates, Barber2016_gaps, Barber2016_costs]= Barber2016(y, ObjectSize, A, M, S, T, NbIters, lambda, theta, TVthresholds, runOnGPU);

% Save the results
filename = ['/home/mory/data/MatlabOneStep/Barber/Barber2016_', num2str(TVthresholds(1)), '_', num2str(TVthresholds(2)), '_' , num2str(TVthresholds(3)), '.mat'];
save(filename, 'Barber2016_iterates', 'Barber2016_gaps', 'Barber2016_costs', '-v7.3');

%% Play iterates

% PlayIterates(Cai2013_iterates, [0 0 0], [0.015 0.015 1.5]);
% PlayIterates(Long2014_iterates, [0.0075 0.0075 0.75], [0.0125 0.0125 1.25]);
% PlayIterates(Weidinger2016_iterates, [0.0075 0.0075 0.75], [0.0125 0.0125 1.25]);
% PlayIterates(Mechlem2017_iterates,  [0.0075 0.0075 0.75], [0.0125 0.0125 1.25]);
% PlayIterates(Barber2016_iterates, [0 0 0], [0.015 0.015 1.5]);