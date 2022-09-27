function [Asol,Fsol,Rsol,hist] = ADJUST(Y,W,k,T,options)
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

 
% get the size of U
[m,n] = size(W);  
[p,c] = size(T);


if nargin < 5
    options = [];
end

epochs  = getoptions(options,'iterMax',1e3); % maximum number of ALS iterations
epIter  = getoptions(options,'innerIter',2); % inner-iterations
qA      = getoptions(options,'qA',1);        % max-value for matrix A
resTol  = getoptions(options,'resTol',1e-6); % tolerance for residual
relTol  = getoptions(options,'relTol',1e-9); % tolerance for residual progress
progTol = getoptions(options,'progTol',1e-6);% tolerance for iterate progress
rho     = getoptions(options,'rho',1e-2);
EQ_A    = getoptions(options,'EQ_A',0);
dynUpd  = getoptions(options,'dynUpd',1);

%% initialization

A0 = rand(n,k);

for i=1:size(A0,1)      % normalize the rows
    A0(i,:) = A0(i,:)/sum(A0(i,:)); 
end

R0 = rand(k,p);
  
%% main loop

% projector for the constraints for A and R
funProjR = @(X) projSimplexR(X,k,p);
funProjA = @(X) projSimplexA(X,qA,EQ_A,n,k);


% optimization options (required for SPG algorithm)
Opt.maxIter     = epIter;     
Opt.memory      = epIter;
Opt.curvilinear = 1;
Opt.progTol     = 1e-9;
Opt.suffDec     = 1e-9;
Opt.testOpt     = 0;
Opt.verbose     = 0;

% initialization
R = R0;
A = A0;
Q = 0*Y;

fprintf('=======================================================================================\n');
fprintf('   ADJUST - Dictionary-based Joint Unmixing of Spectral Tomographic Projections        \n');
fprintf('=======================================================================================\n');
fprintf('image size         : %d x %d \n',sqrt(n),sqrt(n));
fprintf('measurements       : %d x %d (total : %d) \n',size(Y),numel(Y));
fprintf('materials          : %d \n',k);
fprintf('spectral channels  : %d \n',c);
fprintf('Dictionary size    : %d x %d \n',size(T));
fprintf('epochs             : %d (inner_iterations : %d) \n',epochs,Opt.maxIter);
fprintf('acceleration param : %2.2e \n',rho);
fprintf('-----------------------------------------------------------------------\n');
fprintf('%10s %12s %12s %12s %18s \n','Iteration','Residual','Res-prog','Iter-prog','Time');
fprintf('-----------------------------------------------------------------------\n');

normY     = norm(Y,'fro');
time_hist = 0;

% solve multiple times (in case the solution gets stuck)
for i=1:epochs
    
    % solve for R
    tic;
    funObjR = @(X) misfitR(X,Y+rho*Q,W,A,T,k,normY);
    Rp      = R;
    R       = minConf_SPG(funObjR,R(:),funProjR,Opt);
    R       = reshape(R,k,p);
    
    
    % channel matrix
    F       = R*T;
    
    % solve for A
    funObjA = @(X) misfitA(X,F,Y+rho*Q,W,k,normY);
    Ap      = A;
    A       = minConf_SPG(funObjA,A(:),funProjA,Opt);
    A       = reshape(A,n,k);
    
    % update Q
    res  = (Y - W*(A*F));
    Q    = Q + res;
    
    time_iter = toc;
    time_hist = time_hist + time_iter;
    
    % history
    hist.res(i)  = norm(res,'fro')/normY;
    hist.prog(i) = norm(R-Rp,'fro')+norm(A-Ap,'fro');
    hist.rho(i)  = rho;
    if i>1 
        hist.res_prog(i) = abs(hist.res(i) - hist.res(i-1));
    else 
        hist.res_prog(i)=1; 
    end
    
    % show iterate and progress
    fprintf('%10d %12.3e %12.3e %12.3e [%4.2e/%4.2e]\n',i,...
        hist.res(i),hist.res_prog(i),hist.prog(i),time_iter,time_hist);
    
    % record solution
    if i==1
        Asol = A;Fsol = F;Rsol = R;
        hist.resmin = hist.res(i);
    else
       if hist.res(i) <= hist.resmin
          hist.resmin = hist.res(i);
          Asol = A;Fsol = F;Rsol = R;
       end
    end
    
    % stopping criteria
    if (hist.res(i) < resTol)  || (hist.prog(i) < progTol) ...
            || (hist.res_prog(i) < relTol)
        fprintf('%10d %12.3e %12.3e %12.3e\n',i,hist.res(i),...
            hist.res_prog(i),hist.prog(i));
        break;
    end
    
    % dynamic update of rho
    if dynUpd
        if (rho > 5e-2) && (hist.res_prog(i) > 1e-3)
            rho = 1e-2;
            fprintf('Restoring the original rho value \n');
        end
        if (hist.res_prog(i) < 1e-5)
            rho = rho*10;
            fprintf('Increasing rho by factor 10\n');
        end
    end
    
end

fprintf('=======================================================================================\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,g] = misfitR(R,Y,W,A,T,k,nY)
% misfit function for matrix R, where Y = W * A * R * T. 

% get the size of U and reshape matrix R
[m,t] = size(Y);
[p,t] = size(T);
R     = reshape(R,k,p);
F     = R*T;

res = W*(A*F) - Y;                  % residual
f   = 0.5 * norm(res,'fro')^2/nY;   % function value (least-squares misfit)
g   = A' * (W' * (res * T'));       % gradient wrt matrix R (least-squares) 
g   = g(:)/nY;

end

function [f,g] = misfitA(A,F,Y,W,k,nY)
% misfit function handle for matrix A, where Y = W * A * F

% get the size of U and reshape A
[m,n] = size(W);
A     = reshape(A,n,k);

res = W*A*F - Y;                    % residual
f   = 0.5 * norm(res,'fro')^2/nY;   % function value (least-squares misfit)
g   = W' *(res * F');               % gradient wrt A
g   = g(:)/nY;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x] = projSimplexA(x,q,EQ,m,n)
% projector for matrix A (simplex constraints - on rows)

A = reshape(x,m,n);         % reshape matrix R
B = proj_simplex_q(A',q,EQ);% project onto simplex (transpose required!)
x = vec(B');                % transpose and vectorize it

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x] = projSimplexR(x,m,n)

R = reshape(x,m,n);         % reshape matrix R

R = projectTwoSimplexConstr(R);
x = R(:);

end

function [x] = projectTwoSimplexConstr(R)

[m,n]    = size(R);
MAX_ITER = 1e4;

%%% solve using ADMM

x = R;
if sum(x(:)) >= 1
    for i=1:MAX_ITER
        
        x1 = transpose(proj_simplex_q(0.5*(x+R)',1,1));
        x = proj_simplex_q(x1,1,0);
        
        if norm(x1-x,'fro') < 1e-6
            break;
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x] = proj_simplex_q(x,q,EQ)
% simplex projector

% EQ = 0;             % Equality constraints indicator (set to 1 if required)
VECTORIZE = false;  % vectorization indictor ('false' to keep it as matrix)


if any( x(:) < 0 ) || any( sum( x ) > q ) || ( EQ && norm( sum( x ) - q)>1e-8/norm(q) )
    if VECTORIZE
        s     = sort( x(:), 'descend' );
    else
        s     = sort( x, 'descend' );
    end
    if q < eps(s(1))
        error('Input is scaled so large compared to q that accurate computations are difficult');
        % since then cs(1) = s(1) - q is  not even guaranteed to be
        % smaller than q !
    end
    if size(x,2) == 1 || VECTORIZE
        cs    = ( cumsum(s) - q ) ./ ( 1 : numel(s) )';
        ndx   = nnz( s > cs );
        tau   = cs(ndx);
        if ~EQ, tau = max( tau, 0 ); end
        x     = max( x - tau, 0 );
    else
        % vectorized for several columns
        cs    = diag( 1./( 1 : size(s,1) ) )*( cumsum(s) - q );
        ndx   = sum( s > cs );
        ndx   = sub2ind( size(x), ndx, 1:size(x,2)  );
        tau   = cs(ndx);
        if ~EQ, tau = max( tau, 0 ); end
        x     = max( x - repmat(tau,size(x,1),1), 0 );
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = getoptions(options, name, v, mandatory)
% getoptions - retrieve options parameter


if nargin<4
    mandatory = 0;
end

if isfield(options, name)
    v = eval(['options.' name ';']);
elseif mandatory
    error(['You have to provide options.' name '.']);
end

end
