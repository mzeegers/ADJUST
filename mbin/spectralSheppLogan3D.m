function [A] = spectralSheppLogan3D(n)
% SpectralThoraxPhantom3D loads a 3D phantom based on the Shepp-Logan phantom,
% and returns material maps for each material seperately
% (for this phantom the 'background' is included in the returned material
% maps)

    % Load phantom
    I  = phantom3d(n);

    % Get unique grey-levels (number of materials)
    ui = unique(I);
    ui = sort(ui,'descend');
    
    % Number of materials
    m  = length(ui);

    % Initialize material map matrix
    A = zeros(n^3,m-1);

    % Find indices and put them to 1
    for i = 1:m-1
        xId      = find(I==ui(i));
        A(xId,i) = 1;
    end

end
