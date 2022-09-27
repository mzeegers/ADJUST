function toleranceReached = GetFirstIterationMeetingTolerance(iterates, writeEveryNthIterate, tolerance, cropROI)

ObjectSize = sqrt(size(iterates, 1));
unit=ObjectSize/8;

iodinemask = zeros(ObjectSize); iodinemask(2*unit+1+cropROI:3*unit-cropROI, 2*unit+1+cropROI:3*unit-cropROI) = ones(unit-2*cropROI, unit-2*cropROI);
gadomask = zeros(ObjectSize); gadomask(4*unit+1+cropROI:5*unit-cropROI, 5*unit+1+cropROI:6*unit-cropROI) = ones(unit-2*cropROI, unit-2*cropROI);
watermask = zeros(ObjectSize); watermask(unit+1+cropROI:7*unit-cropROI, unit+1+cropROI:7*unit-cropROI) = ones(6*unit-2*cropROI, 6*unit-2*cropROI);

maskediodine = squeeze(iterates(:,1,:)) .* iodinemask(:); meanIodine = squeeze(sum(maskediodine, 1)) / sum(iodinemask(:));
maskedgado = squeeze(iterates(:,2,:)) .* gadomask(:); meanGado = squeeze(sum(maskedgado, 1)) / sum(gadomask(:));
maskedwater = squeeze(iterates(:,3,:)) .* watermask(:); meanWater = squeeze(sum(maskedwater, 1)) / sum(watermask(:));

targetiodine = 0.01;
targetgado = 0.01;
targetwater = 1;

% Compute when all the materials have reached the target concentration
% up to the requested tolerance
inBounds = ((abs(meanIodine - 0.01) < (targetiodine * tolerance)) ...
         .* (abs(meanGado   - 0.01) < (targetgado   * tolerance)) ...
         .* (abs(meanWater  - 1)    < (targetwater  * tolerance)));

% Get the first iteration from which the concentrations are in the bounds
% and never leave them anymore
toleranceReached = numel(inBounds) - find(fliplr(1 - inBounds), 1);

if (numel(toleranceReached) == 0)
    toleranceReached = 0;
else
    toleranceReached = toleranceReached * writeEveryNthIterate;
end

end