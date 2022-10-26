function [Asol, Fsol, hist] = cJoint(Y, W, k, options)
% cJoint - classical JOINT reconstruction and unmixing algorithm for
% solving the spectral tomographic inverse problem
%
%               Y = W * A * F
%
% by posing it as a following optimization problem
%
%  minimize_{A,F}  || W * A * F - Y ||_F^2
%  subject to      A >= 0
%                  F >= 0
%
%  W is a tomography projection matrix (size: m x n)
%  Y is the set projections over all frequencies (size: m x c)
%  m is the number of projections
%  c is the number of frequencies
%  n is the number of voxels/pixels (image size)
%  A has a size n x k (with k << n)
%  F has a size k x c
%  ||X||_F = \sqrt{ \sum_{i=1}^{n} (x_i)^2 } denotes Frobenius norm
%
% Input:
%   Y - spectral tomographic projections (size: m x c)
%   W - tomography forward projection operator (size: m x n)
%   k - number of materials in the object
%   options - a strcture containing following options
%       epochs: number of outer iterations (default: 1e3)
%       epIter: number of inner iterations (default: 2)
%       resTol: residual tolerance for optimality condition (default: 1e-6)
%       progTol:iterate progress tolerance (default: 1e-6)
%       rho   : acceleration parameter for ADJUST (default: 1e-2)
%
% Output:
%   Asol - (optimal) material map from cJoint (size: n x k)
%   Fsol - (optimal) spectral map from cJoint (size: k x c)
%   hist - history structure containing following
%       res : residual for each iterate (size: t x 1)
%       prog: progress of variables at each iterate (size: t x 1)
%       rho : acceleration parameter for each iterate (size: t x 1)
%       here t is the number of cJoint iterations
%
% Authors:
%   Ajinkya Kadu,
%       Centrum Wiskunde & Informatica, Amsterdam (aak@cwi.nl)
%   Mathé Zeegers, 
%       Centrum Wiskunde & Informatica, Amsterdam (M.T.Zeegers@cwi.nl)

if nargin<4, options = []; end

epochs  = getoptions(options, 'iterMax' ,  1e3);   % maximum number of iterations
epIter  = getoptions(options, 'innerIter', 3);     % inner-iterations
resTol  = getoptions(options, 'resTol',    1e-6);  % tolerance for residual
progTol = getoptions(options, 'progTol',   1e-6);  % tolerance for progress
rho     = getoptions(options, 'rho',       0.01);  % acceleration parameter
showIter= getoptions(options, 'showIter',  1);     % display iter
use_cuda= getoptions(options, 'use_cuda',  gpuDeviceCount > 0); % use GPUs

%% initialization

% get the size of U
[m, n] = size(W);
[m, c] = size(Y);

% randomized spatial maps
A0 = rand(n, k);

% initial spectral maps
F0 = rand(k, c);

if use_cuda
    A0 = gpuArray(single(A0));
    F0 = gpuArray(single(F0));
    Y  = gpuArray(single(Y));
end

%% main loop

% projector for the constraints for A and R
funProjF = @(X) projF(X);
funProjA = @(X) projA(X);

% optimization options (required for SPG algorithm)
Opt.maxIter     = epIter;
Opt.memory      = epIter;
Opt.curvilinear = 1;
Opt.progTol     = 1e-9;
Opt.suffDec     = 1e-9;
Opt.testOpt     = 0;
Opt.verbose     = 0;

% initialization
F = F0;
A = A0;
Q = 0*Y;

fprintf('==========================================================\n');
fprintf('                cJoint                                    \n');
fprintf('==========================================================\n');
fprintf('image size         : %d x %d \n',sqrt(n),sqrt(n));
fprintf('measurements       : %d x %d (total : %d) \n',size(Y),numel(Y));
fprintf('materials          : %d \n',k);
fprintf('spectral channels  : %d \n',c);
fprintf('epochs             : %d (inner_iterations : %d) \n',epochs,Opt.maxIter);
fprintf('-----------------------------------------------------------------------\n');
fprintf('%10s %12s %12s %12s %18s \n','Iteration','Residual','Res-prog','Iter-prog','Time');
fprintf('-----------------------------------------------------------------------\n');

normY     = norm(Y,'fro');
time_hist = 0;

% solve multiple times (in case the solution gets stuck)
for i=1:epochs

    % solve for F
    tic;
    funObjF = @(X) misfitF(X, Y + rho*Q , W, A, k, normY, use_cuda);
    Fp      = F;
    F       = minConf_SPG(funObjF, F(:), funProjF, Opt);
    F       = reshape(F, k, c);

    % solve for A
    funObjA = @(X) misfitA(X, F, Y + rho*Q, W, k, normY, use_cuda);
    Ap      = A;
    A       = minConf_SPG(funObjA, A(:), funProjA, Opt);
    A       = reshape(A, n, k);

    % update Q
    if use_cuda
        res = Y - gpuArray(W * gather(A)) * F;
    else
        res  = Y - (W * A) * F;
    end
    Q    = Q + res;

    % register time
    time_iter = toc;
    time_hist = time_hist + time_iter;

    % history
    hist.res(i)  = norm(res, 'fro') / normY;
    hist.prog(i) = norm(F - Fp, 'fro') + norm(A - Ap, 'fro');
    hist.rho(i)  = rho;
    if i>1
        hist.res_prog(i) = hist.res(i) - hist.res(i-1);
    else
        hist.res_prog(i)=1; 
    end

    % show iterate and progress
    if ((showIter == 1) && (rem(i,100)==0)) || (showIter > 1)
    fprintf('%10d %12.3e %12.3e %12.3e [%4.2e/%4.2e]\n',i,...
        hist.res(i),hist.res_prog(i),hist.prog(i),time_iter,time_hist);
    end

    % record solution with least residual (since the convegence is non-monotonic)
    if i==1
        Asol = A;
        Fsol = F;
        hist.resmin = hist.res(i);
    else
        if hist.res(i) < hist.resmin
            Asol = A;
            Fsol = F;
            hist.resmin = hist.res(i);
        end
    end

    % stopping criteria
    if (hist.res(i) < resTol)  || (hist.prog(i) < progTol)
        fprintf('%10d %15.5e %15.5e\n', i, hist.res(i), hist.prog(i));
        break;
    end

end

if use_cuda
    Asol = gather(Asol);
    Fsol = gather(Fsol);
end

fprintf('==========================================================\n\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,g] = misfitF(F, Y, W, A, k, nY, use_cuda)
% misfit function for matrix F, where Y = W * A * F.
% provide functional value and gradient at given point

% get the size of U and reshape matrix R
[m, c] = size(Y);
F      = reshape(F, k, c);

if use_cuda
    WA  = gpuArray(W * gather(A));
else
    WA  = W * A;
end

res = WA * F - Y;                   % residual
f   = 0.5 * norm(res,'fro')^2 / nY; % function value (least-squares misfit)
g   = WA' * res;                    % gradient wrt matrix F (least-squares)
g   = g(:) / nY;

end


function [f,g] = misfitA(A, F, Y, W, k, nY, use_cuda)
% misfit function handle for matrix A
% provide functional value and gradient at given point

% get the size of U and reshape A
[m, n] = size(W);
A      = reshape(A, n, k);

if use_cuda
    WA  = gpuArray(W * gather(A));
    res = WA * F - Y;                   % residual
    f   = 0.5 * norm(res,'fro')^2 /nY;  % function value (least-squares misfit)
    g   = W' * gather(res * F');        % gradient wrt A
    g   = gpuArray( g(:) )/ nY;
else
    res = (W*A) * F - Y;                % residual
    f   = 0.5 * norm(res,'fro')^2 /nY;  % function value (least-squares misfit)
    g   = W' * (res * F');              % gradient wrt A
    g   = g(:) / nY;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x] = projA(x)
% projector for matrix A (bounds constraints)

x(x < 0) = 0;
x(x > 1) = 1;

end

function [x] = projF(x)
% project x onto non-negative orthant

x(x < 0) = 0;

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
