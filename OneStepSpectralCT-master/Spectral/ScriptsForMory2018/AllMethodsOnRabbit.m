%% Load all data

% Import photon counts
pc = read_mhd('/home/mory/data/Hamburg/data_2018_Feb_13.9714.1/pcounts.mhd');
measuredCounts = pc.data;
[pix_per_proj, rows, projs, bins] = size(measuredCounts);
measuredCounts = reshape(measuredCounts, [pix_per_proj * projs, bins]); % Fails if rows ~= 1

% Import detector response and spectrum model
det = read_mhd('/home/mory/data/Hamburg/realScannerModel/detector_response.mhd');
detectorResponse = det.data';
[detectedEnergies, incidentEnergies] = size(detectorResponse);

sp = read_mhd('/home/mory/data/Hamburg/data_2018_Feb_13.9714.1/spectrumForMatlab.mhd');
incidentSpectrum = squeeze(sp.data).';
% IMPORTANT: rtkprojectionmatrix does not accept projections with non-identity transform. 
% Since PRH projections come in with negative spacing on x, they had to be clitkAffineTransform-ed
% and re-ordered. The spectrum therefore needs to be reordered too
% So this is only for data whose projection matrix is exported via
% rtkprojectionmatrix
incidentSpectrum = fliplr(incidentSpectrum);

mat = read_mhd('/home/mory/data/Hamburg/realScannerModel/material_attenuations.mhd');
M = mat.data';

S_unbinned = detectorResponse .* reshape(incidentSpectrum, [1 incidentEnergies pix_per_proj]);

% Define thresholds and get S binned
thresholds = [30; 51; 64; 72; 85];
S = zeros(numel(thresholds), incidentEnergies, pix_per_proj);
thresholds = cat(1, thresholds, size(S_unbinned, 1));

shiftThresholds=false;
if (shiftThresholds)
    for i=1:size(S, 1)
        S(i,:,:) = sum(S_unbinned(thresholds(i)-3:thresholds(i+1)-2,:,:), 1);
    end
else
    for i=1:size(S, 1)
        S(i,:,:) = S_unbinned(thresholds(i), :,:) / 2 + sum(S_unbinned(thresholds(i)+1:thresholds(i+1)-1,:,:), 1) + S_unbinned(thresholds(i+1), :,:) / 2;
    end
end

% % Import sparse geometry matrix, transform it into what Matlab needs
% fileID = fopen('/home/mory/data/Hamburg/data_2018_Feb_13.9714.1/sparseProjectionMatrix.bin');
% A = fread(fileID,'double');
% A = reshape(A, [3,numel(A)/3]);
% A(1,:) = A(1,:)+1;
% A(2,:) = A(2,:)+1;
% A = spconvert(A');
% fclose(fileID);
% 
% % Save the matrix once transformed, then the next times just reload it
% save('/home/mory/data/Hamburg/data_2018_Feb_13.9714.1/sparseProjectionMatrix.mat', 'A', '-v7.3');

% Load the matrix
load('/home/mory/data/Hamburg/data_2018_Feb_13.9714.1/sparseProjectionMatrix.mat', 'A');

% Set parameters that are identical for all methods
ObjectSize = 380;
runOnGPU = true;
delta_huber = [0.001, 0.001, 0.01];

% y = measuredCounts;
% save('/home/mory/data/MatlabOneStep/preps_rabbit.mat', 'y', 'A', 'M', 'S', '-v7.3');

%% Reconstruct with SQS methods

% No change of basis, use T=identity
T = eye(size(M, 2));

% Run Mechlem2017
NbIters = 200;
NbSubsets = 4; % More than that speeds up computation, but makes it unstable, and ultimately, divergent
% regulWeights = [30000, 30000, 10]; % Huber is scaled by delta, so where there is a low delta, there should be a high lambda to compensate
regulWeights = [10000, 10000, 10]; % Huber is scaled by delta, so where there is a low delta, there should be a high lambda to compensate
[Mechlem2017_iterates, Mechlem2017_costs ]= Mechlem2017(measuredCounts, ObjectSize, pix_per_proj, A, M, S, T, NbIters, regulWeights, NbSubsets, delta_huber, runOnGPU);

% Save the results
filename = ['/home/mory/data/MatlabOneStep/Mechlem/Rabbit_Mechlem2017_', num2str(regulWeights(1)), '_', num2str(regulWeights(2)), '_' , num2str(regulWeights(3)), '.mat'];
save(filename, 'Mechlem2017_iterates', 'Mechlem2017_costs', '-v7.3');

% Run Weidinger2016
NbIters = 5000;
regulWeights = [3000, 3000, 1]; %Perfect
[Weidinger2016_iterates, Weidinger2016_costs ]= Weidinger2016(measuredCounts, ObjectSize, A, M, S, T, NbIters, regulWeights, runOnGPU);

% Save the results
filename = ['/home/mory/data/MatlabOneStep/Weidinger/Rabbit_Weidinger2016_', num2str(regulWeights(1)), '_', num2str(regulWeights(2)), '_' , num2str(regulWeights(3)), '.mat'];
save(filename, 'Weidinger2016_iterates', 'Weidinger2016_costs', '-v7.3');

% Run Long2014
NbIters = 5000;
NbSubsets = 20; % For ordered subsets
% regulWeights = [30000, 30000, 300]; % Way too much regul
% regulWeights = [3000, 3000, 10]; % Not enough regul
% regulWeights = [10000, 10000, 30]; % Good, but too much noise on contrast
% images
regulWeights = [30000, 30000, 30];
[Long2014_iterates, Long2014_costs ]= Long2014(measuredCounts, ObjectSize, pix_per_proj, A, M, S, T, NbIters, regulWeights, NbSubsets, delta_huber, runOnGPU);

% Save the results
filename = ['/home/mory/data/MatlabOneStep/Long/Rabbit_Long2014_', num2str(regulWeights(1)), '_', num2str(regulWeights(2)), '_' , num2str(regulWeights(3)), '.mat'];
save(filename, 'Long2014_iterates', 'Long2014_costs', '-v7.3');

%% Barber2016

% Compute normalized mu-preconditioning
[M_norm, T_norm] = mu_preconditioning(M, 'normalize', S);

% Run Barber2016
NbIters = 10000;
lambda = 0.0001; % Controls convergence
theta = 0.5; % Between 0 and 1, see Chambolle Pock
TVthresholds = [200,200,20000];
[Barber2016_iterates, Barber2016_gaps, Barber2016_costs]= Barber2016(measuredCounts, ObjectSize, A, M_norm, S, T_norm, NbIters, lambda, theta, TVthresholds, runOnGPU);

% Save the results
filename = ['/home/mory/data/MatlabOneStep/Barber/Rabbit_Barber2016_', num2str(TVthresholds(1)), '_', num2str(TVthresholds(2)), '_' , num2str(TVthresholds(3)), '.mat'];
save(filename, 'Barber2016_iterates', 'Barber2016_gaps', 'Barber2016_costs', '-v7.3');

%% Reconstruct with Cai

% Compute Fessler's mu-preconditioning
normalized_S = S ./ (sum(S, 2) + eps);
[M_fessler, T_fessler] = mu_preconditioning(M, 'fessler', normalized_S);

% Normalize the photon counts by the expected counts if no object
cellS = cell(size(S,3), 1);
for pix=1:size(S,3)
    cellS{pix} = sparse(squeeze(S(:,:,pix)));
end
diagS = blkdiag(cellS{:});
NbPixels = size(A, 1);
[NbBins, NbEnergies, NbPixelsPerProj] = size(S);
NbProjs = NbPixels / NbPixelsPerProj;
NbMats = size(M,2);
zeroprojs = A * zeros(ObjectSize^2, NbMats);
attenuationFactors = exp( - M * zeroprojs.');
reshaped_att = reshape(attenuationFactors, [NbEnergies * NbPixelsPerProj NbProjs]);
forw = diagS * reshaped_att;
forw = reshape(forw, [NbBins NbPixels]).';
y = measuredCounts ./ (forw + eps);
k_d = 0.003;

% Run Cai2013
NbIters = 5000;
useNesterov = false;
% regulWeights = [1000, 1000, 1]; Not enough regul
regulWeights = [10000, 10000, 10];
[Cai2013_iterates, Cai2013_costs ]= Cai2013(y, ObjectSize, A, M_fessler, normalized_S, k_d, NbIters, T_fessler, regulWeights, delta_huber, runOnGPU, useNesterov);

% Save the results
filename = ['/home/mory/data/MatlabOneStep/Cai/Rabbit_Cai2013_', 'fessler', '_', num2str(regulWeights(1)), '_', num2str(regulWeights(2)), '_' , num2str(regulWeights(3)), '.mat'];
save(filename, 'Cai2013_iterates', 'Cai2013_costs', '-v7.3');