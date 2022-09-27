function [sep] = separationCurves(input)
%SEPARATIONINDEX returns a measure of how well the materials of input are
%separated
%   Compute the correlation coefficients matrix, and averages the
%   non-diagonal elements

[~, ~, niter] = size(input);
sep = zeros(1,niter);
for it = 1:niter
    sep(1, it) = separationIndex(input(:,:,it));
end

figure; plot(sep);
end

