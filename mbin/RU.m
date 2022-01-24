function [A,F] = RU(Y,W,k)
%RU Reconstruction-then-Unmixing method for spectral tomography
% 
% We solve the problem 
%                       W * A * F = Y 
% by two-step approach
%   (Reconstruction)    W * U = Y 
%   (Unmixing)          A * F = Y
%
% Input:
%   Y - Spectral tomographic projections (size: m x c)
%   W - Tomographic projection operator (size: m x n)
%   k - Number of materials (note: k << min(m,c))
%
% Here, m is number of projections per spectral bin, c is the number of
% spectral bins, n is the size of image (i.e., total number of
% voxels/pixels)
%
% Output:
%   A - spatial distribution of materials (size: n x k)
%   F - spectral distribution of materials (size: k x c)
%
% Authors:
%   Ajinkya Kadu,
%       Centrum Wiskunde & Informatica, Amsterdam (aak@cwi.nl)
%   Mathé Zeegers, 
%       Centrum Wiskunde & Informatica, Amsterdam (M.T.Zeegers@cwi.nl)

% get sizes
[m,c] = size(Y);

%% main loop
fprintf('==========================================================\n');
fprintf('     Reconstruction - then - Unmixing  \n');
fprintf('==========================================================\n');

%%% reconstruction with CG and Tikhonov regularization
regT = 1e-3;        % Tikhonov reg parameter 
Wty  = W'*Y;
WtW  = @(x) W'*(W*x) + regT*x;

tol  = 1e-6;
iter = 20;
for i=1:c       % run for each channel
   [Up(:,i),~] = pcg(WtW,Wty(:,i),tol,iter);
end

% make spectral volume non-negative
Up(Up <0) = 0;

%%% unmixing - Non-negative least squares (Mult followed by ALS)
opt = statset('MaxIter',10,'Display','final');
[A10,F10] = nnmf(Up,k,'Options',opt,'Algorithm','mult','Replicates',10);

opt = statset('MaxIter',1000,'Display','final');
[A,F] = nnmf(Up,k,'W0',A10,'H0',F10,'Options',opt,'Algorithm','als','Replicates',10);

fprintf('==========================================================\n');

end
