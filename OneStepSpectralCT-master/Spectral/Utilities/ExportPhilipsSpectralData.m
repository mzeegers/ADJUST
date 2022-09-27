function out = ExportPhilipsSpectralData(input, outputDir)

% %MaxEnergy = max([max(input.spectrum.energy) max(input.detector.response.energy) max(input.material.iodine.E) max(input.material.gadolinium.E) max(input.material.water.E)])
MaxEnergy = 150;
% cylindricalDetectorRadius = input.source_rotation_axis + input.detector_rotation_axis;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Export photon counts
% temp = input.photon_count_data;
% NDims = numel(size(temp.data)) - 1;
% 
% % Permute to put the energy as the first dimension
% temp.data = permute(temp.data, [NDims+1 1:NDims]);
% temp.size = size(temp.data).';
% 
% % Permute, to match the usual RTK orientation of projections
% temp.data = permute(temp.data, [1, 3, 2, 4]);
% t = temp.size(2);
% temp.size(2) = temp.size(3);
% temp.size(3) = t;
% 
% % Compute the spacing from fan and cone beam angle lists
% % The fan and cone beam angles point towards the centers of the pixels, 
% % therefore the "-1"
% verticalSpacing = (input.spectrum.tags.SpectrumConeAngleMaximum - input.spectrum.tags.SpectrumConeAngleMinimum) / (input.spectrum.tags.NumberSpectrumConeAngles - 1) * pi / 180 * cylindricalDetectorRadius;
% horizontalSpacing = (input.spectrum.tags.SpectrumFanAngleMaximum - input.spectrum.tags.SpectrumFanAngleMinimum) / (input.spectrum.tags.NumberSpectrumFanAngleValues - 1) * pi / 180 * cylindricalDetectorRadius;
% temp.spacing = ones(NDims, 1);
% temp.spacing(1) = horizontalSpacing;
% temp.spacing(2) = verticalSpacing;
% 
% % Compute the origin, assuming that the detector is centered. Offsets
% % will be accounted for in the geometry description
% origin = zeros(NDims, 1);
% origin(1) = (input.spectrum.tags.NumberSpectrumFanAngleValues - 1) * horizontalSpacing / 2;
% origin(2) = - (numel(input.spectrum.cone_angle) - 1) * verticalSpacing / 2;
% orientation = eye(NDims);
% orientation(1,1) = -1;
% 
% temp.origin = origin;
% temp.orientation = orientation;
% % export_vector_mhd(strcat(outputDir, '/photon_count.mhd'), temp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export PRH decomposed projections
NDims = 3;
siz = zeros(NDims, 1);
origin = zeros(NDims, 1);
orientation = eye(NDims);
orientation(1,1) = -1;
horizontalSpacing = 0.4902;
verticalSpacing = 0.5;
siz(1) = input.proj_gadolinium_medianfiltered.tags.NumberDetectorColumns;
siz(2) = input.proj_gadolinium_medianfiltered.tags.NumberDetectorRows;
siz(3) = input.proj_gadolinium_medianfiltered.tags.NumberTotalViews;

PRH_proj = ImageType(siz,origin,[horizontalSpacing; verticalSpacing; 1],orientation);
% PRH_proj.data = permute(input.proj.gd.data, [2 1 3]);
% write_mhd(strcat(outputDir, '/PRH_proj_gd.mhd'), PRH_proj);
PRH_proj.data = permute(input.proj_gadolinium_medianfiltered.data, [2 1 3]);
write_mhd(strcat(outputDir, '/PRH_proj_gd_median_filtered.mhd'), PRH_proj);

% PRH_proj.data = permute(input.proj.iodine.data, [2 1 3]);
% write_mhd(strcat(outputDir, '/PRH_proj_iodine.mhd'), PRH_proj);
PRH_proj.data = permute(input.proj_iodine_medianfiltered.data, [2 1 3]);
write_mhd(strcat(outputDir, '/PRH_proj_iodine_median_filtered.mhd'), PRH_proj);

% PRH_proj.data = permute(input.proj.water.data, [2 1 3]);
% write_mhd(strcat(outputDir, '/PRH_proj_water.mhd'), PRH_proj);
PRH_proj.data = permute(input.proj_water_medianfiltered.data, [2 1 3]);
write_mhd(strcat(outputDir, '/PRH_proj_water_median_filtered.mhd'), PRH_proj);
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export incident spectrum
temp = input.spectrum;
NDims = numel(size(temp.data)) - 1;

% Partly take the 'energy' field into account by shifting the table
tempData = zeros(MaxEnergy, size(temp.data, 2), size(temp.data, 3));
tempData(temp.energy(1):temp.energy(1)+size(temp.data, 1)-1, :, :) = temp.data;
temp.data = tempData;

% Take into account mAs, solid angle and exposure time
solidAngle = verticalSpacing * horizontalSpacing / (input.source_rotation_axis + input.detector_rotation_axis)^2;
temp.data = temp.data * input.mA * solidAngle * input.rotation_time/input.number_views_per_turn;
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
export_vector_mhd(strcat(outputDir, '/spectrum.mhd'), temp);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Export detector response as an image
% temp = input.detector.response;
% 
% % Partly take the 'energy' field into account by shifting the table
% tempData = zeros(size(temp.data, 1), MaxEnergy);
% tempData(:, temp.energy(1):temp.energy(1)+size(temp.data, 2)-1) = temp.data;
% temp.data = tempData;
% 
% % Detector response matrix
% drm = temp.data;
% 
% temp.data = cast(temp.data.', 'single');
% temp.size = size(temp.data).';
% temp.origin = zeros(numel(temp.size), 1);
% temp.spacing = temp.origin + 1;
% temp.orientation = eye(numel(temp.size));
% write_mhd(strcat(outputDir, '/detector_response.mhd'), temp);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Export material attenuations as an image
% temp = input.material;
% 
% NbMaterials = numel(fieldnames(temp));
% 
% % Partly take the 'energy' field into account by shifting the table
% tempData = zeros(NbMaterials, MaxEnergy);
% 
% materialNames = fieldnames(temp);
% for i=1:NbMaterials
%    materialName = materialNames{i};
%    materialStruct = getfield(temp, materialName);
%    tempData(i, materialStruct.E(1):materialStruct.E(1)+size(materialStruct.E, 2)-1) = materialStruct.mu;
% end
% temp.data = tempData;
% temp.data = cast(temp.data, 'single');
% temp.size = size(temp.data).';
% temp.origin = zeros(numel(temp.size), 1);
% temp.spacing = temp.origin + 1;
% temp.orientation = eye(numel(temp.size));
% write_mhd(strcat(outputDir, '/material_attenuations.mhd'), temp);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Export reconstructed volumes
% temp.data = permute(input.img.iodine.data, [2 3 1]);
% temp.size = size(temp.data);
% temp.spacing = [0.3906 0.3906 0.3906];
% temp.origin = -(temp.size - 1) .* temp.spacing / 2;
% temp.orientation = eye(numel(temp.size));
% write_mhd(strcat(outputDir, '/img_iodine.mhd'), temp);
% 
% temp.data = permute(input.img.gd.data, [2 3 1]);
% write_mhd(strcat(outputDir, '/img_gd.mhd'), temp);
% 
% temp.data = permute(input.img.water.data, [2 3 1]);
% write_mhd(strcat(outputDir, '/img_water.mhd'), temp);

% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Heuristic for decomposition starting point estimation
% 
% % Compute the input spectrum as seen by the detector
% wspdatamatrix = reshape(wspdata, [size(wspdata, 1), size(wspdata, 2) * size(wspdata, 3)]);
% isp_asbd = drm * wspdatamatrix;
% isp_asbd = reshape(isp_asbd, [size(isp_asbd, 1), size(wspdata, 2), size(wspdata, 3)]);
% 
% % Bin it
% binned_isp_asbd = zeros(numel(thresholds), size(isp_asbd, 2), size(isp_asbd, 3));
% for bin = 1:numel(thresholds)-1
%     binned_isp_asbd(bin, :,:) = sum(isp_asbd(thresholds(bin):thresholds(bin+1)-1, :,:), 1);
% end
% binned_isp_asbd(numel(thresholds), :,:) = sum(isp_asbd(thresholds(numel(thresholds)):size(isp_asbd, 1), :,:), 1);
% 
% % Log transform the photon counts
% logs = -log(input.photon_count_data.data ./ permute(repmat(binned_isp_asbd, [1 1 1 2400]), [2 3 4 1]));
% 
% % Compute the mean attenuation of each material in each bin
% mean_att_binned = zeros(size(tempData, 1), numel(thresholds));
% for bin = 1:numel(thresholds)-1
%     mean_att_binned(:,bin) = mean(tempData(:,thresholds(bin):thresholds(bin+1)-1), 2);
% end
% mean_att_binned(:,numel(thresholds)) = mean(tempData(:,thresholds(numel(thresholds)):size(tempData, 2)), 2);
% 
% % Estimate the length of each material traversed as if all attenuation had
% % to be explained by that material
% estimatedLengths = repmat(logs, [1 1 1 1 NbMaterials]) ./ repmat(reshape(mean_att_binned.', [1 1 1 size(mean_att_binned, 2), size(mean_att_binned, 1)]), [size(logs, 1), size(logs, 2), size(logs, 3), 1, 1]);
% startingPoint = squeeze(min(estimatedLengths, [], 4));
% % startingPoint = squeeze(estimatedLengths(:,:,:,3,:));
% 
% % Export it like the photon count data
% temp = input.photon_count_data;
% temp.data = startingPoint;
% NDims = numel(size(temp.data)) - 1;
% 
% % Permute to put the energy as the first dimension
% temp.data = permute(temp.data, [NDims+1 1:NDims]);
% temp.size = size(temp.data).';
% 
% % Permute, to match the usual RTK orientation of projections
% temp.data = permute(temp.data, [1, 3, 2, 4]);
% t = temp.size(2);
% temp.size(2) = temp.size(3);
% temp.size(3) = t;
% temp.spacing = ones(NDims, 1);
% temp.spacing(1) = horizontalSpacing;
% temp.spacing(2) = verticalSpacing;
% 
% % Compute the origin, assuming that the detector is centered. Offsets
% % will be accounted for in the geometry description
% origin = zeros(NDims, 1);
% origin(1) = (input.spectrum.tags.NumberSpectrumFanAngleValues - 1) * horizontalSpacing / 2;
% origin(2) = - (numel(input.spectrum.cone_angle) - 1) * verticalSpacing / 2;
% orientation = eye(NDims);
% orientation(1,1) = -1;
% 
% temp.origin = origin;
% temp.orientation = orientation;
% export_vector_mhd(strcat(outputDir, '/starting_point.mhd'), temp);


out = 'SUCCESS';