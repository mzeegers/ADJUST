function out = ExportPhilipsDetectorResponseFromBrownEtAl(input, lowerDetectorPart, upperDetectorPart, outputDir)

MaxEnergy = 150;
cylindricalDetectorRadius = input.source_rotation_axis + input.detector_rotation_axis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export photon counts
temp = input.photon_count_data;
NDims = numel(size(temp.data)) - 1;

% Permute to put the energy as the first dimension
temp.data = permute(temp.data, [NDims+1 1:NDims]);
temp.size = size(temp.data).';

% Permute, to match the usual RTK orientation of projections
temp.data = permute(temp.data, [1, 3, 2, 4]);
t = temp.size(2);
temp.size(2) = temp.size(3);
temp.size(3) = t;

% Compute the spacing from fan and cone beam angle lists
% The fan and cone beam angles point towards the centers of the pixels, 
% therefore the "-1"
verticalSpacing = (input.spectrum.tags.SpectrumConeAngleMaximum - input.spectrum.tags.SpectrumConeAngleMinimum) / (input.spectrum.tags.NumberSpectrumConeAngles - 1) * pi / 180 * cylindricalDetectorRadius;
horizontalSpacing = (input.spectrum.tags.SpectrumFanAngleMaximum - input.spectrum.tags.SpectrumFanAngleMinimum) / (input.spectrum.tags.NumberSpectrumFanAngleValues - 1) * pi / 180 * cylindricalDetectorRadius;
temp.spacing = ones(NDims, 1);
temp.spacing(1) = horizontalSpacing;
temp.spacing(2) = verticalSpacing;

% Compute the origin, assuming that the detector is centered. Offsets
% will be accounted for in the geometry description
origin = zeros(NDims, 1);
origin(1) = (input.spectrum.tags.NumberSpectrumFanAngleValues - 1) * horizontalSpacing / 2;
origin(2) = - (numel(input.spectrum.cone_angle) - 1) * verticalSpacing / 2;
orientation = eye(NDims);
orientation(1,1) = -1;

temp.origin = origin;
temp.orientation = orientation;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export the product between input spectrum and detector response, for the
% lower part of the detector (high energies)
temp = input.spectrum;
NDims = numel(size(temp.data)) - 1;
temp.data = repmat(lowerDetectorPart.', [1 size(temp.data, 2) size(temp.data, 3)]);

% Take into account mAs, solid angle and exposure time
solidAngle = verticalSpacing * horizontalSpacing / (input.source_rotation_axis + input.detector_rotation_axis)^2;
temp.data = temp.data * input.mA * solidAngle * input.rotation_time/input.number_views;
temp.data = cast(temp.data, 'single');
temp.size = size(temp.data).';

% Permute, to match the usual RTK orientation of projections
temp.data = permute(temp.data, [1, 3, 2]);
t = temp.size(2);
temp.size(2) = temp.size(3);
temp.size(3) = t;

temp.origin = zeros(NDims, 1);
temp.spacing = temp.origin + 1;
temp.orientation = eye(NDims);
export_vector_mhd(strcat(outputDir, '/lowerDetectorResponse.mhd'), temp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export the product between input spectrum and detector response, for the
% upper part of the detector (low energies)
temp = input.spectrum;
NDims = numel(size(temp.data)) - 1;
temp.data = repmat(upperDetectorPart.', [1 size(temp.data, 2) size(temp.data, 3)]);

% Take into account mAs, solid angle and exposure time
solidAngle = verticalSpacing * horizontalSpacing / (input.source_rotation_axis + input.detector_rotation_axis)^2;
temp.data = temp.data * input.mA * solidAngle * input.rotation_time/input.number_views;
temp.data = cast(temp.data, 'single');
temp.size = size(temp.data).';

% Permute, to match the usual RTK orientation of projections
temp.data = permute(temp.data, [1, 3, 2]);
t = temp.size(2);
temp.size(2) = temp.size(3);
temp.size(3) = t;

temp.origin = zeros(NDims, 1);
temp.spacing = temp.origin + 1;
temp.orientation = eye(NDims);
export_vector_mhd(strcat(outputDir, '/upperDetectorResponse.mhd'), temp);


out = 'SUCCESS';