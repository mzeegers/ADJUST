function [A] = spectralThoraxPhantom()
% SpectralThoraxPhantom loads the Thorax phantom,
% and returns material maps for each material seperately 

    % Load phantom
    I  = load('ThoraxPhantom512Slice255Mod.mat').J;

    % Get unique grey-levels (number of materials)
    ui = unique(I);
    ui = sort(ui,'descend');
    m  = length(ui);

    % Initialize material map matrix (exclude background)
    A = zeros(size(I,1)*size(I,2),m-1);

    % Find indices and put them to 1
    for i = 1:m-1
        xId      = find(I==ui(i));
        A(xId,i) = 1;
    end

end
