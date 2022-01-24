function [A,X] = spectralPhantomCirclesMixed(N,k)
% SpectralPhantomCircles computes and loads the Mixed Disk phantom (solid disks and mixed disks),
% and returns material maps for each material seperately
% Note: Use maximum k=5 (or k=6 but there are some circles with cuts) and minimum k=2

    % Initialize material map matrix (excluding background)
    A = zeros(k,N, N);

    % Fixed radius for each disk
    radius = round(N/15);
        
    % Put in each solid disk
    for i=0:k-1
        
        % Find center of the i-th disk (use a reduced radius to put these on the inner circle)
        centery = round((3/8)*(N/2)* cos(2*pi*i/k) + (N/2));
        centerx = round((3/8)*(N/2)* sin(2*pi*i/k) + (N/2));
        
        %Now create the disk with the given radius
        x = -1*radius:radius;
        y = -1*radius:radius;
        [xx yy] = meshgrid(x,y);
        u = zeros(size(xx));
        u((xx.^2+yy.^2)<(radius)^2)=1;
        
        % Insert into the phantom
        A(i+1,centery-radius:centery+radius,centerx-radius:centerx+radius) = u;
        
    end

    %Now add in the mixed disks (on the outer circle
    %Determine the number of disks K
    K = nchoosek(k,2);
    mixone = 1;
    mixtwo = 2;

    for i=0:K-1
        
        %Find center of the i-th disk
        %Use a reduced radius
        centery = round((3/4)*(N/2)* cos(2*pi*i/K) + (N/2));
        centerx = round((3/4)*(N/2)* sin(2*pi*i/K) + (N/2));

        %Now create the disk with the given radius
        x = -1*radius:radius;
        y = -1*radius:radius;
        [xx yy] = meshgrid(x,y);
        u = zeros(size(xx));
        u((xx.^2+yy.^2)<(radius)^2)=1;
        
        %Add the mixture of two different disks
        A(mixone,centery-radius:centery+radius,centerx-radius:centerx+radius) = 0.5*u;
        A(mixtwo,centery-radius:centery+radius,centerx-radius:centerx+radius) = 0.5*u;
        
        if(mixtwo == k)
            mixone = mixone + 1;
            mixtwo = mixone;
        end
        mixtwo = mixtwo + 1;
        
    end

    % Reshape the array to return flattened material maps
    A = reshape(A, k, N*N);
    A = A.';

end
