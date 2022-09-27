function [toplot, iterations, finalmean, finalstd] = getMaterialConcentration(material, iterates, writeEveryNthIterate, cropROI)
%GETMATERIALCONCENTRATION Outputs material's mean concentration over iterations
% Make the ROI smaller by cropROI pixels than the reference square containing each material, to ignore border effects

ObjectSize = sqrt(size(iterates, 1));
unit=ObjectSize/8;

if (strcmp(material,'iodine'))
    mask = zeros(ObjectSize); mask(2*unit+1+cropROI:3*unit-cropROI, 2*unit+1+cropROI:3*unit-cropROI) = ones(unit-2*cropROI, unit-2*cropROI);
    maskedIterates = squeeze(iterates(:,1,:)) .* mask(:);
end
if (strcmp(material,'gadolinium'))
    mask = zeros(ObjectSize); mask(4*unit+1+cropROI:5*unit-cropROI, 5*unit+1+cropROI:6*unit-cropROI) = ones(unit-2*cropROI, unit-2*cropROI);
    maskedIterates = squeeze(iterates(:,2,:)) .* mask(:);
end
if (strcmp(material,'water'))
    mask = zeros(ObjectSize); mask(unit+1+cropROI:7*unit-cropROI, unit+1+cropROI:7*unit-cropROI) = ones(6*unit-2*cropROI, 6*unit-2*cropROI);
    maskedIterates = squeeze(iterates(:,3,:)) .* mask(:);
end

toplot = squeeze(sum(maskedIterates, 1)) / sum(mask(:));
finalmean = toplot(end);
tmp = maskedIterates(:,end);
finalstd = std(tmp(mask(:)>0));

iterations = (1:numel(toplot)) * writeEveryNthIterate;
end