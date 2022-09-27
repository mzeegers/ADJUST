function [A,b,x,theta,p,d] = odftomo(N,theta,p,d)
%ODFTOMO Creates a 2D tomography test problem with a smooth odf image
%
%   [A b x theta p d] = odftomo(N)
%   [A b x theta p d] = odftomo(N,theta)
%   [A b x theta p d] = odftomo(N,theta,p)
%   [A b x theta p d] = odftomo(N,theta,p,d)
%
% This function creates a 2D tomography test problem with an N-times-N
% domain which represents a smooth odf image of a porous material, using
% p parallel rays for each angle in the vector theta.
%
% Input: 
%   N           Scalar denoting the number of discretization intervals in 
%               each dimesion, such that the domain consists of N^2 cells.
%   theta       Vector containing the angles in degrees. Default: theta = 
%               0:5:175.
%   p           Number of parallel rays for each angle. Default: p = N.
%   d           Scalar denoting the distance from the first ray to the last.
%               Default: d = N.
%
% Output:
%   A           Coefficient matrix with N^2 columns and nA*p rows, 
%               where nA is the number of angles, i.e., length(theta).
%   b           Vector containing the rhs of the test problem.
%   x           Vector containing the exact solution, with elements
%               between 0 and 1.
%   theta       Vector containing the used angles in degrees.
%   p           The number of used rays for each angle.
%   d           The distance between the first and the last ray.
% 
% See also: binarytomo, fanbeamtomo, paralleltomo, seismictomo.

% Per Christian Hansen, June 8, 2012, DTU Informatics.

% Default value of d.
if nargin < 5 || isempty(seed), seed = 17; end

% Default value of d.
if nargin < 4 || isempty(d), d = N; end

% Default value of the number of rays.
if nargin < 3 || isempty(p), p = round(2*N); end

% Default value of the angles theta.
if nargin < 2 || isempty(theta), theta = 0:5:175; end

% Generate the image.
[I,J] = meshgrid(1:N);
sigma = 0.25*N;
c = 0.6*N;
X = exp( - (I-c).^2/(1.2*sigma)^2 - (J-c).^2/sigma^2) + ...
    0.5*exp( - (I-0.5*N).^2/(1.2*sigma)^2 - (J-0.3*N).^2/sigma^2);
x = X(:);

A = paralleltomo(N,theta,p,d);
b = A*x;