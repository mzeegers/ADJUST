function [Asol, Fsol, Rsol, hist] = ADJUST(Y, W, k, T, options)
% ADJUST - A Dictionary-based Joint Unmixing and reconstruction for
% Spectral Tomography. Our approach solves the problem 
%
%               Y = W * A * F
%
% by incorporating a spectral dictionary T and representing the spectral
% signatures using F = R * T, where R is a dictionary coefficient matrix.
% The resulting problem is
%  
%  minimize_{A,R}  || W * A * R * T - Y ||_F^2 
%  subject to      A \in row_simplex
%                  R \in {row_simplex, column_simplex} 
%
%  W is a tomography projection matrix (size: m x n)
%  Y is the set projections over all frequencies (size: m x c)
%  m is the number of projections
%  c is the number of frequencies
%  n is the number of voxels/pixels (image size)
%  T is a dictionary matrix consists of materials with channel frequency
%       information (size: p x c)
%  A has a size n x k (with k << n)
%  R has a size k x p 
%  ||X||_F = \sqrt{ \sum_{i=1}^{n} (x_i)^2 } denotes Frobenius norm
%  row_simplex means that the sum of rows should be less than equal to 1
%       and all elements must be greater than 0.
%  column_simplex means the sum of elements along column should be less
%       than or equal to 1 and all elements must be > 0.
%
% Input:
%   Y - spectral tomographic projections (size: m x c)
%   W - tomography forward projection operator (size: m x n)
%   k - number of materials in the object
%   T - spectral dictionary (size: p x c)
%   options - a strcture containing following options
%       epochs: number of outer iterations (default: 1e3)
%       epIter: number of inner iterations (default: 2)
%       qA    : maximum value for matrix A (default: 1)
%       resTol: residual tolerance for optimality condition (default: 1e-6)
%       relTol: relative residual tolerance (default: 1e-9)
%       progTol:iterate progress tolerance (default: 1e-6)
%       rho   : acceleration parameter for ADJUST (default: 1e-2)
%       EQ_A  : (strict) equality constraint, on/off (default: 0)
%       dynUpd: experimental parameter to dynamically update rho values
%               (default: 1)
%
% Output:
%   Asol - (optimal) material map from ADJUST (size: n x k)
%   Fsol - (optimal) spectral map from ADJUST (size: k x c)
%   Rsol - (optimal) spectral dictionary coefficients (size: k x p)
%   hist - history structure containing following
%       res : residual for each iterate (size: t x 1)
%       prog: progress of variables at each iterate (size: t x 1)
%       rho : acceleration parameter for each iterate (size: t x 1)
%       here t is the number of ADJUST iterations
%
% Authors:
%   Ajinkya Kadu,
%       Centrum Wiskunde & Informatica, Amsterdam (aak@cwi.nl)
%   Mathé Zeegers, 
%       Centrum Wiskunde & Informatica, Amsterdam (M.T.Zeegers@cwi.nl)

if nargin < 5, options = []; end

epochs  = getoptions(options, 'iterMax',  1e3);  % maximum number of ALS iterations
epIter  = getoptions(options, 'innerIter',5);    % inner-iterations
qA      = getoptions(options, 'qA',       1);    % max-value for matrix A
resTol  = getoptions(options, 'resTol',   1e-9); % tolerance for residual
relTol  = getoptions(options, 'relTol',   1e-9); % tolerance for residual progress
progTol = getoptions(options, 'progTol',  1e-9); % tolerance for iterate progress
rho     = getoptions(options, 'rho',      1e-2); % acceleration parameter
EQ_A    = getoptions(options, 'EQ_A',     0);    % equality constraints on A (1: strong, 0: weak)
dynUpd  = getoptions(options, 'dynUpd',   1);    % dynamic update of acceleration parameter (experimental)
showIter= getoptions(options, 'showIter', 1);    % show iterates (1: yes, 0: no)
use_cuda= getoptions(options, 'use_cuda', gpuDeviceCount > 0); % Use GPU arrays         

%% initialization

% get the size of U
[m, n] = size(W);  
[p, c] = size(T);

A0 = rand(n, k);
R0 = rand(k, p);

for i = 1:size(A0,1)      % normalize the rows by their sum
    A0(i,:) = A0(i,:) / sum(A0(i,:)); 
end

%%% convert to GPU arrays
if use_cuda
    A0 = gpuArray(single(A0));
    R0 = gpuArray(single(R0));
    Y  = gpuArray(single(Y));
    T  = gpuArray(single(T));
end
  
%% main loop

% projector for the constraints for A and R
funProjR = @(X) projSimplexR(X, k, p);
funProjA = @(X) projSimplexA(X, qA, EQ_A, n, k);

% optimization options (required for SPG algorithm)
Opt.maxIter     = epIter;     
Opt.memory      = epIter;
Opt.curvilinear = 1;
Opt.progTol     = 1e-12;
Opt.suffDec     = 1e-12;
Opt.testOpt     = 0;
Opt.verbose     = 0;

% initialization
R = R0;
A = A0;
Q = 0*Y;

fprintf('==========================================================\n');
fprintf('    ADJUST - Dictionary-based Joint Unmixing              \n');
fprintf('             of Spectral Tomographic Projections           \n');
fprintf('==========================================================\n');
fprintf('image size         : %d x %d \n', sqrt(n), sqrt(n));
fprintf('measurements       : %d x %d (total : %d) \n', size(Y), numel(Y));
fprintf('materials          : %d \n', k);
fprintf('spectral channels  : %d \n', c);
fprintf('Dictionary size    : %d x %d \n', size(T));
fprintf('epochs             : %d (inner_iterations : %d) \n', epochs, Opt.maxIter);
fprintf('acceleration param : %2.2e \n', rho);
fprintf('-----------------------------------------------------------------------\n');
fprintf('%10s %12s %12s %12s %18s \n','Iteration','Residual','Res-prog','Iter-prog','Time');
fprintf('-----------------------------------------------------------------------\n');

normY     = norm(Y,'fro');
time_hist = 0;

% solve multiple times (in case the solution gets stuck)
for i = 1:epochs
    
    % solve for R
    tic;
    funObjR = @(X) misfitR(X, Y + rho*Q, W, A, T, k, normY, use_cuda);
    Rp      = R;
    R       = minConf_SPG(funObjR, R(:), funProjR, Opt);
    R       = reshape(R, k, p);
    
    % channel matrix
    F       = R*T;
    
    % solve for A
    funObjA = @(X) misfitA(X, F, Y + rho*Q, W, k, normY, use_cuda);
    Ap      = A;
    A       = minConf_SPG(funObjA, A(:), funProjA, Opt);
    A       = reshape(A, n, k);
    
    % update Q
    if use_cuda
        res  = Y - gpuArray(W * gather(A)) * F;
    else
        res  = Y - (W * A) * F;
    end
    Q    = Q + res;
    
    % register time
    time_iter = toc;
    time_hist = time_hist + time_iter;
    
    % history
    hist.res(i)  = norm(res, 'fro') / normY;
    hist.prog(i) = norm(R - Rp, 'fro') + norm(A - Ap, 'fro');
    hist.rho(i)  = rho;
    if i > 1 
        hist.res_prog(i) = abs(hist.res(i) - hist.res(i-1));
    else 
        hist.res_prog(i) = 1; 
    end
    
    % show iterate and progress (every 100 iterations)
    if ((showIter == 1) && (rem(i,100)==0)) || (showIter > 1)
        fprintf('%10d %12.3e %12.3e %12.3e [%4.2e/%4.2e]\n',i,...
            hist.res(i),hist.res_prog(i),hist.prog(i),time_iter,time_hist);
    end
    
    % record solution
    if i==1
        Asol = A;
        Fsol = F;
        Rsol = R;
        hist.resmin = hist.res(i);
    else
       if hist.res(i) <= hist.resmin
          Asol = A;
          Fsol = F;
          Rsol = R;
          hist.resmin = hist.res(i);
       end
    end
    
    % stopping criteria
    if (hist.res(i) < resTol)  || (hist.prog(i) < progTol) ...
            || (hist.res_prog(i) < relTol)
        fprintf('%10d %12.3e %12.3e %12.3e\n',i,hist.res(i),...
            hist.res_prog(i),hist.prog(i));
        break;
    end
    
    % dynamic update of rho (this is an experimental feature)
    if dynUpd
        if (rho > 5e-2) && (hist.res_prog(i) > 1e-3)
            rho = 1e-2;
            if showIter > 1
                fprintf('Restoring the original rho value \n'); 
            end
        end
        if (hist.res_prog(i) < 1e-5)
            rho = rho * 10;
            if showIter > 1, fprintf('Increasing rho by factor 10\n'); end
        end
    end
    
end

if use_cuda
    Asol = gather(Asol);
    Fsol = gather(Fsol);
    Rsol = gather(Rsol);
end

fprintf('==========================================================\n\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, g] = misfitR(R, Y, W, A, T, k, nY, use_cuda)
% misfit function for matrix R, where Y = W * A * R * T. 
% provides functional value and gradient at particular point

% get the size of U and reshape matrix R
[m, t] = size(Y);
[p, t] = size(T);
R      = reshape(R, k, p);
F      = R * T;

if use_cuda
    WA  = gpuArray(W * gather(A));
else
    WA  = W * A;
end

res = WA * F - Y;                    % residual
f   = 0.5 * norm(res, 'fro')^2 / nY; % function value (least-squares misfit)
g   = WA' * (res * T');              % gradient wrt matrix R (least-squares)
g   = g(:) / nY;

end

function [f, g] = misfitA(A, F, Y, W, k, nY, use_cuda)
% misfit function handle for matrix A, where Y = W * A * F
% provides functional value and gradient at particular point

% get the size of U and reshape A
[m, n]= size(W);
A     = reshape(A, n, k);

if use_cuda
    WA  = gpuArray(W * gather(A));
    res = WA * F - Y;                     % residual
    f   = 0.5 * norm(res, 'fro')^2 / nY;  % function value (least-squares misfit)
    g   = gpuArray(W' * gather(res * F'));% gradient wrt A
    g   = g(:) / nY;
else
    res = (W * A) * F - Y;              % residual
    f   = 0.5 * norm(res, 'fro')^2 / nY;% function value (least-squares misfit)
    g   = W' * (res * F');              % gradient wrt A
    g   = g(:) / nY;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x] = projSimplexA(x, q, EQ, m, n)
% projector for matrix A (simplex constraints - on rows)

A = reshape(x, m, n);          % reshape matrix R
B = proj_simplex_q(A', q, EQ); % project onto simplex (transpose required!)
x = vec(B');                   % transpose and vectorize it

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x] = projSimplexR(x, m, n)
% projector for matrix R (simplex constraints on both rows and columns)

R = reshape(x, m, n);         % reshape matrix R
R = projectTwoSimplexConstr(R);
x = R(:);

end

function [x] = projectTwoSimplexConstr(R)

[m,n]    = size(R);
MAX_ITER = 1e4;

%%% solve using alternating projections (similar to Dijkstra's algorithm)
x = R;

if sum(x(:)) >= 1
    for i = 1 : MAX_ITER
        
        % row projection
        % x1 = proj_simplex_q(0.5 * (x + R)', 1, 0);
        x1 = proj_simplex_q(x', 1, 1);
        x1 = x1';
        
        % column projection
        x = proj_simplex_q(x1, 1, 0);
        
        % stopping criterion
        if norm(x1 - x, 'fro') < 1e-6, break; end
        
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x] = proj_simplex_q(x, q, EQ)
% simplex projector
% EQ: Equality constraints indicator (set to 1 if required)


if any( x(:) < 0 ) || any( sum(x) > q ) || ( EQ && norm(sum(x) - q) > 1e-6 / norm(q) )
    s     = sort( x, 'descend' );
    if q < eps(s(1))
        error('Input is scaled so large compared to q that accurate computations are difficult');
        % since then cs(1) = s(1) - q is  not even guaranteed to be
        % smaller than q !
    end
    % vectorized for several columns
    cs    = diag( 1./( 1 : size(s, 1) ) ) * ( cumsum(s) - q );
    ndx   = sum( s > cs );
    ndx   = sub2ind( size(x), ndx, 1:size(x,2) );
    tau   = cs(ndx);
    if ~EQ, tau = max( tau, 0 ); end
    x     = max( x - repmat(tau, size(x, 1), 1), 0 );
end

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
