function [ diffs ] = SQSdiffWithNeighbors( input )
%SQSDIFFWITHNEIGHBORS Computes the difference with surrounding pixels
% Assumes the input image is 3D (2D, channels), and returns 
% a 4D image (2D, channels, neighboring pixel). Each of the (:,:,:,?)
% slices contains the difference between the center pixel and one of its
% neighbors

N = 3;
[height, width, channels] = size(input);
diffs = zeros(height, width, channels, N^2, 'like', input);

padpattern = [(N-1)/2 (N-1)/2 0];
paddedinput = padarray(input, padpattern, 'both', 'replicate');
count=0;
for i=1:N
    for j=1:N
        count = count+1;
        diffs(:,:,:,count) = input - paddedinput(i:i+height-1, j:j+width-1, :);
    end
end



end