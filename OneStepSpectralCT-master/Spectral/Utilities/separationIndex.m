function [sep] = separationIndex(input)
%SEPARATIONINDEX returns a measure of how well the materials of input are
%separated
%   Compute the correlation coefficients matrix, and averages the
%   non-diagonal elements

[~, nmat] = size(input);
stdev = std(input) + eps;
std_prods = stdev' * stdev;
correlationCoeffs = cov(input) ./ std_prods;
nondiag = correlationCoeffs - eye(nmat);
sep = sum(abs(nondiag(:))) / (numel(nondiag) - nmat);
end

