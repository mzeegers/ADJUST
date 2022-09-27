function out = ExportDucrosSpectralData(fwd, rec, outputDir)

MaxEnergy = 120;
MinEnergy = 12;
NDims = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export photon counts

% % Permute, to match the usual RTK orientation of projections
% temp.data = permute(temp.data, [1, 3, 2, 4]);
% t = temp.size(2);
% temp.size(2) = temp.size(3);
% temp.size(3) = t;

temp.spacing = ones(NDims, 1);
temp.origin =  zeros(NDims, 1);
temp.orientation = eye(NDims);
temp.data = reshape(fwd.S0, [4, 256, 256]);
temp.size = size(temp.data).';
export_vector_mhd(strcat(outputDir, '/noiseless_photon_count.mhd'), temp);

temp.data = reshape(fwd.S, [4, 256, 256]);
export_vector_mhd(strcat(outputDir, '/noisy_photon_count.mhd'), temp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export incident spectrum

% Pad
spectrum = zeros(MaxEnergy, 1);
spectrum(MinEnergy:MaxEnergy) = fwd.N0;

% Make the spectrum an image
spectrum = repmat(squeeze(spectrum), [1 256 256]);
temp.data = spectrum;
temp.size = size(temp.data).';
export_vector_mhd(strcat(outputDir, '/spectrum.mhd'), temp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export detector response as an image

% Pad
response = zeros(4,MaxEnergy);
response(:,MinEnergy:MaxEnergy) = fwd.D;

% temp.data = S;
temp.data = response.';
temp.size = size(temp.data).';
write_mhd(strcat(outputDir, '/detector_response.mhd'), temp);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export material attenuations as an image

%Pad
attenuations = zeros(MaxEnergy, 3);
attenuations(MinEnergy:MaxEnergy, :) = rec.T;

temp.data = attenuations.';
temp.size = size(temp.data).';
write_mhd(strcat(outputDir, '/material_attenuations.mhd'), temp);

out = 'SUCCESS';