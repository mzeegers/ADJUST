% This script simulates preps through a simple phantom, then attempts to 
% reconstruct the phantom with the Cai2013 method, internally using a basis of
% synthetic materials to speed up convergence

% Set the size of the phantom
ObjectSize = 256;
MaxIters = 200;
NbPixelsPerProj = round(ObjectSize * sqrt(2));

% Simulate noisy preps
noise = true;
getBeamTransmissionRatios = true;

%[y, k_d, A, M, S, T, ~] = SimulatePreps(ObjectSize, noise, true, 'fessler');

% attenuations
load('../ScannerModel/incidentSpectrum.mat', 'incidentSpectrum');
load('../ScannerModel/materialAttenuations.mat', 'materialAttenuations');
load('../ScannerModel/detectorResponse.mat', 'detectorResponse');

% Get them under matrix form, consistent with Cai's paper
materialAttenuations = materialAttenuations(25:124,:);
M = materialAttenuations;
%S_unbinned = detectorResponse .* incidentSpectrum.';
S_unbinned = incidentSpectrum.';

% Define thresholds and get S binned
%S = cat(1, sum(S_unbinned(30:50, :), 1), ...
%           sum(S_unbinned(51:61, :), 1), ...
%           sum(S_unbinned(62:71, :), 1), ...
%           sum(S_unbinned(72:82, :), 1), ...
%           sum(S_unbinned(83:end, :), 1));

S = S_unbinned(:);
S = S(25:124);
S = diag(S);


% Generate a simple phantom
unit = ObjectSize/8;
water=zeros(ObjectSize);
water(unit+1:7*unit, unit+1:7*unit) = ones(6*unit, 6*unit);
iodine=zeros(ObjectSize);
iodine(2*unit+1:3*unit, 2*unit+1:3*unit) = ones(unit, unit) * 0.01 / 4.933;
gadolinium=zeros(ObjectSize);
gadolinium(4*unit+1:5*unit, 5*unit+1:6*unit) = ones(unit, unit) * 0.01 / 7.9;

% writeImageResult(iodine, 'iodineGT.png', 0, 0.01);
% writeImageResult(gadolinium, 'gadoliniumGT.png', 0, 0.01);
% writeImageResult(water, 'waterGT.png', 0, 1);

% Compute the number of projections required
% to keep a pixels-to-voxels ratio of 4
NbProj = ceil(ObjectSize * 4 / sqrt(2));
ForwardTheta = [0:NbProj-1] * 180 / NbProj;

% Generate forward projection matrix
A = paralleltomo(ObjectSize,ForwardTheta);

% Compute forward projections
x = cat(2, iodine(:), gadolinium(:), water(:));
%
tic

object = x;
noisy = noise;
normalize = getBeamTransmissionRatios;

% Compute forward projections
projs = A * object;

% Compute photon counts as seen by the detector
attenuationFactors = exp(- M * projs.');
counts = S * attenuationFactors;
counts_without_attenuation = S * ones(size(attenuationFactors));

if (noisy)
	% Add noise to the photon counts
	noisy_counts = poissrnd(counts);
else
	noisy_counts = counts;
end

if (normalize)
	% Divide by the counts in the absence of object
	y = noisy_counts ./ counts_without_attenuation;
else
	y = noisy_counts;
end

y_without_attenuation = poissrnd(counts_without_attenuation) ./ counts_without_attenuation;
k_d = var(y_without_attenuation(:)) / mean(y_without_attenuation(:));

y = y.';


if (getBeamTransmissionRatios)
    % Normalize the joint spectra matrix S
    noAttenuation = S * ones(size(M,1),1);
    S = S ./ repmat(noAttenuation, [1, size(S, 2)]);
end
toc

S = repmat(S, [1 1 NbPixelsPerProj]); % Same spectrum, repeated for all pixels. Just to test for pixel-dependent spectrum feature

% Set regularization parameters
lambda = [1, 1, 1];

% Perform reconstruction
runOnGPU = false;
[Weidinger2016_iterates, costs ]= Weidinger2016(y, ObjectSize, A, M, S, T, NbIters, lambda, runOnGPU);

% Multiply the results by each material's density, in order to get something in g/ml
Weidinger2016_iterates(:,1,:)=Weidinger2016_iterates(:,1,:)*4.933; %Iodine
Weidinger2016_iterates(:,2,:)=Weidinger2016_iterates(:,2,:)*7.9; %Gadolinium
% Nothing to do for water, since it has density 1

% Show the resulting sequence of iterates
PlayIterates(Weidinger2016_iterates, [0 0 0], [0.015 0.015 1.5])

% Plot the resulting costs, in loglog scale
loglog(costs - min(costs(:)))
