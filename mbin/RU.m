function [A,F] = RU(Y, W, k, options)
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

if nargin < 4, options = []; end

regT    = getoptions(options, 'regParam',   1e-3); % Tikhonov regularization parameter
pcgIter = getoptions(options, 'pcgIter',    20);   % PCG max iterations
pcgTol  = getoptions(options, 'pcgTol',     1e-6); % PCG tolerance
multIter= getoptions(options, 'multIter',   100);  % NMF mult algorithm max iterations
alsIter = getoptions(options, 'alsIter',    30);   % NMF als algorithm max iterations
reps    = getoptions(options, 'replicates', 10);   % number of replicates for NMF
usePar  = getoptions(options, 'useParallel',license('test','Distrib_Computing_Toolbox'));

%% main loop

fprintf('==========================================================\n');
fprintf('     Reconstruction - then - Unmixing  (RU)               \n');
fprintf('==========================================================\n');

%%% reconstruction with CG and Tikhonov regularization
Wty  = W'*Y;
WtW  = @(x) W' * (W * x) + regT * x;


fprintf('running reconstruction (per channel)... \n');
Up = zeros(size(Wty));

for i = 1:c       % run for each channel

    [Up(:,i), ~] = pcg(WtW, Wty(:,i), pcgTol, pcgIter);
   
    % print progress
    if mod(i, floor(0.2*c)) == 1, fprintf('%d ... ', i); ...
    elseif i == c, fprintf('%d \n', i); end
end

% make spectral volume non-negative
Up(Up < 0) = 0;

%%% unmixing - Non-negative least squares (Mult followed by ALS)
opt = statset('MaxIter', multIter, 'Display', 'final', ... 
                'UseParallel', usePar);
[A10, F10] = nnmf(Up, k, 'Options', opt, 'Algorithm', 'mult', ...
                'Replicates', reps);

opt = statset('MaxIter', alsIter, 'Display', 'final', ...
                'UseParallel', usePar);
[A, F] = nnmf(Up, k, 'W0', A10, 'H0', F10, 'Options', opt, ...
                'Algorithm', 'als', 'Replicates', reps);

fprintf('==========================================================\n\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = getoptions(options, name, v, mandatory)
% getoptions - retrieve options parameter

if nargin < 4, mandatory = 0; end

if isfield(options, name)
    v = eval(['options.' name ';']);
elseif mandatory
    error(['You have to provide options.' name '.']);
end

end
