function [A, F] = UR(Y, W, k, options)
%UR Unmixing-then-Reconstruction for spectral tomography
% 
% We solve the problem 
%                       W * A * F = Y 
% by two-step approach
%   (Unmixing)          P * F = Y
%   (Reconstruction)    W * A = P 
%
% Input:
%   Y - Spectral tomographic projections (size: m x c)
%   W - Tomographic projection operator (size: m x n)
%   k - Number of materials (note: k << min(m,c))
%
% Here, m is number of projections per spectral bin, c is the number of
% spectral bins, n is the size of image
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
[m,n] = size(W);

if nargin < 4, options = []; end

regT    = getoptions(options, 'regParam',   1e-3); % Tikhonov regularization parameter
pcgIter = getoptions(options, 'pcgIter',    100);  % PCG max iterations
pcgTol  = getoptions(options, 'pcgTol',     1e-6); % PCG tolerance
multIter= getoptions(options, 'multIter',   100);  % NMF mult algorithm max iterations
alsIter = getoptions(options, 'alsIter',    30);   % NMF als algorithm max iterations
reps    = getoptions(options, 'replicates', 10);   % number of replicates for NMF
usePar  = getoptions(options, 'useParallel',license('test','Distrib_Computing_Toolbox'));

%% main loop

fprintf('==========================================================\n');
fprintf('     Unmixing - then - Reconstruction (UR)                \n');
fprintf('==========================================================\n');

%%% unmixing - Non-negative least squares (Mult followed by ALS)
opt = statset('MaxIter', multIter, 'Display', 'final', ...
             'UseParallel', usePar);
[A0, F0] = nnmf(Y, k, 'Options', opt, 'Algorithm', 'mult', ...
             'replicates', reps);

opt = statset('MaxIter', alsIter, 'Display', 'final', ...
              'UseParallel', usePar);
[P, F] = nnmf(Y, k, 'W0', A0, 'H0', F0, 'Options', opt, ...
              'algorithm', 'als', 'Replicates', reps);


%%% reconstruction with CG and Tikhonov regularization
Wtp  = W' * P;      % backprojection
WtW  = @(x) W' *( W * x) + regT * x;  % regularized forward operator

fprintf('running reconstruction (per material)... \n');

A    = zeros(n, k);

for i = 1:k           % run for each material

    [A(:,i), ~] = pcg(WtW, Wtp(:,i), pcgTol, pcgIter);

    % print progress
    if i < k, fprintf('%d ... ', i); else, fprintf('%d \n', i); end
end

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

