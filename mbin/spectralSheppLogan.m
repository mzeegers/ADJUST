function [A] = spectralSheppLogan(n)
% SpectralThoraxPhantom loads the Shepp-Logan phantom,
% and returns material maps for each material seperately
% (for this phantom the 'background' is included in the returned material
% maps)

    % Load phantom
    I = phantom(n);
    
    % Vectorize the array
    I = abs(I(:));

    % Get unique grey-levels (number of materials)
    ui = unique(I);
    ui = sort(ui, 'descend');

    % Simplify the phantom to 5 materials (including background)
    I(I == ui(2)) = ui(3);
    I(I == ui(5)) = ui(3);

    ui = unique(I);
    ui = sort(ui, 'descend');

    % Number of materials
    m = length(ui);

    % Initialize material map matrix (including background)
    A = zeros(n^2, m);

    % Find indices and put them to 1
    for i = 1:m
        xId = find(I == ui(i));
        A(xId, i) = 1;
    end

end
