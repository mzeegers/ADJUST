function [A, X] = spectralPhantomCircles(N, k)
% SpectralPhantomCircles computes and loads the Disk phantom,
% and returns material maps for each material seperately

    % Initialize material map matrix (excluding background)
    A = zeros(k, N, N);

    % Fixed radius for each disk
    radius = round(N / 15);
    
    % Put in each disk
    for i = 0:k-1
        
        % Find center of the i-th disk
        centery = round((3/4) * (N/2) * cos(2*pi*i/k) + (N/2));
        centerx = round((3/4) * (N/2) * sin(2*pi*i/k) + (N/2));

        % Now create the disk with the given radius
        radius = round(radius);
        x = -1*radius:radius;
        y = -1*radius:radius;
        [xx, yy] = meshgrid(x, y);
        u = zeros(size(xx));
        u(xx.^2 + yy.^2 < radius^2) = 1;
           
        % Insert into the phantom
        A(i+1, centery-radius:centery+radius, ...
            centerx-radius:centerx+radius) = u;
        
    end
   

    % Reshape the array to return flattened material maps
    A = reshape(A, k, N * N);
    A = A.';

end
