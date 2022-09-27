function [] = InvertCovarianceMatrix( input, output )
%INVERTCOVARIANCEMATRIX Summary of this function goes here
%   Detailed explanation goes here

toInvert = read_mhd(input);
A = toInvert.data;
D = zeros(size(A));

% Invert sub-matrices one by one, and check their eigenvalues
for x=1:size(A,1)
    for y=1:size(A,2)
        for z=1:size(A,3)
            mat = reshape(A(x,y,z,:), [sqrt(size(A,4)) sqrt(size(A,4))]).';
            E = eig(mat);
            if (sum(E <= 0) > 0)
                mat
                x
                y
                z
            end
            invmat = inv(mat);
            D(x,y,z,:) = invmat(:);
        end
    end
end
inverted = toInvert;
inverted.data = D;
write_mhd(output, inverted);

