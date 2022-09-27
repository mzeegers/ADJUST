% Import photon counts
pc = read_mhd('/home/mory/data/tests/rtksparseprojectionmatrixexport/photonCounts.mhd');
measuredCounts = pc.data;
[pix_per_proj, rows, projs, bins] = size(measuredCounts);
measuredCounts = reshape(measuredCounts, [pix_per_proj * projs, bins]); % Fails if rows ~= 1

% Import detector response and spectrum model
det = read_mhd('/home/mory/data/Hamburg/realScannerModel/detector_response.mhd');
detectorResponse = det.data.';
[detectedEnergies, incidentEnergies] = size(detectorResponse);

sp = read_mhd('/home/mory/data/Hamburg/2017_Jan_17.7905.1_materialDecomp/spectrumForMatlab.mhd');
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
thresholds = [30; 51; 62; 72; 83];
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

% No change of basis, use T=identity
T = eye(size(M, 2));
       
% % Import sparse geometry matrix
% fileID = fopen('/home/mory/data/tests/rtksparseprojectionmatrixexport/projectionMatrix.bin');
% A = fread(fileID,'double');
% A = reshape(A, [3,numel(A)/3]);
% A(1,:) = A(1,:)+1;
% A(2,:) = A(2,:)+1;
% A = spconvert(A');
% fclose(fileID);
load('/home/mory/data/tests/rtksparseprojectionmatrixexport/forwardProjectionMatrix.mat', 'A');

%% Set reconstruction parameters with Mechlem
ObjectSize = 380;
NbIters = 50;
runOnGPU = true;
NbSubsets = 4; % More than that speeds up computation, but makes it unstable, and ultimately, divergent
regulWeights = [30000, 30000, 30]; % Huber is scaled by delta, so where there is a low delta, there should be a high lambda to compensate
delta_huber = [0.001, 0.001, 0.01];

% Run Mechlem2017
[Mechlem2017_iterates, Mechlem2017_costs ]= Mechlem2017(measuredCounts, ObjectSize, pix_per_proj, A, M, S, T, NbIters, regulWeights, NbSubsets, delta_huber, runOnGPU);

% Multiply by the density
Mechlem2017_iterates(:,1,:)=Mechlem2017_iterates(:,1,:)*4.933; % iodine
Mechlem2017_iterates(:,2,:)=Mechlem2017_iterates(:,2,:)*7.9; % gadolinium
Mechlem2017_iterates(:,3,:)=Mechlem2017_iterates(:,3,:)*1; % water

% View results
PlayIteratesTransposed(Mechlem2017_iterates, [0 0 0], [0.015 0.015 1.5]);
result = Mechlem2017_iterates(:,:,NbIters);
result = reshape(result, [ObjectSize ObjectSize 3]);
mhdstructure = read_mhd('/home/mory/data/tests/rtksparseprojectionmatrixexport/refvol.mhd');
mhdstructure.data = 1000*result(:,:,1); write_mhd('/home/mory/data/MatlabOneStep/Mechlem/phantom_iodine.mhd', mhdstructure);
mhdstructure.data = 1000*result(:,:,2); write_mhd('/home/mory/data/MatlabOneStep/Mechlem/phantom_gadolinium.mhd', mhdstructure);
mhdstructure.data = 1000*result(:,:,3); write_mhd('/home/mory/data/MatlabOneStep/Mechlem/phantom_water.mhd', mhdstructure);

% Compare with PRH-reconstructions
load('/home/mory/data/Hamburg/2017_Jan_17.7905.1_materialDecomp/PhilipsReconstruction/2017_Jan17.7905.1_materialImages.mat', 'materialImages');
PRH_result = cat(3, materialImages.iodine(:,:,4), materialImages.gadolinium(:,:,4), materialImages.water(:,:,4));
PlayIteratesTransposed(reshape(PRH_result, [512*512 3 1]), [0 0 0], [0.015 0.015 1.5]);
mhdstructure = read_mhd('/home/mory/data/MatlabOneStep/Scanner/refvol.mhd');
mhdstructure.data = 1000*materialImages.iodine(:,:,4); write_mhd('/home/mory/data/MatlabOneStep/Scanner/phantom_iodine.mhd', mhdstructure);
mhdstructure.data = 1000*materialImages.gadolinium(:,:,4); write_mhd('/home/mory/data/MatlabOneStep/Scanner/phantom_gadolinium.mhd', mhdstructure);
mhdstructure.data = 1000*materialImages.water(:,:,4); write_mhd('/home/mory/data/MatlabOneStep/Scanner/phantom_water.mhd', mhdstructure);
