function [Asol,Fsol,Rsol,hist] = ADJUST(Y,W,k,T,A0,R0,options)
%  
%  minimize_{R,T}  || W * A * R * T - Y ||^2 
%  subject to      0 <= A_i <= 1
%                  R \in row_simplex
%
%  W is a tomography projection matrix (size: m x n)
%  Y are the projections for each frequency (size: m x k)
%  m are the number of projections
%  n are the number of frequencies
%  T is a dictionary matrix consists of materials with channel frequency
%       information (size: p x k)
%  A has a size n x k (with k << n)
%  R has a size k x p 
%  m_Cone - indicator for monotone cone constraints
%  m_eq - indicator for (strict) equality in simplex constraints
 
% get the size of U
[m,n] = size(W);  
[p,t] = size(T);


if nargin<7
    options = [];
end

epochs  = getoptions(options,'iterMax',1e3); % maximum number of ALS iterations
epIter  = getoptions(options,'innerIter',2); % inner-iterations
qA      = getoptions(options,'qA',1);        % max-value for matrix A
resTol  = getoptions(options,'resTol',1e-6); % tolerance for residual
relTol  = getoptions(options,'relTol',1e-9); % tolerance for residual progress
progTol = getoptions(options,'progTol',1e-6);% tolerance for iterate progress
rho     = getoptions(options,'rho',0);
simFlag = getoptions(options,'simFlag',1);
EQ_A    = getoptions(options,'EQ_A',0);
lambda  = getoptions(options,'lambda',0);
printFig= getoptions(options,'printFig',0);
yWeight = getoptions(options,'yWeight',0);

%%
A0 = rand(n,k);

for i=1:size(A0,1)      % normalize the rows
    A0(i,:) = A0(i,:)/sum(A0(i,:)); 
end

R0 = rand(k,p);
  
%% solve for F (A is solved internally)

if yWeight
    mY= mean(Y,1);
    Omega = spdiags(1./mY',0,size(Y,2),size(Y,2));
    Y = Y*Omega;
    T = T*Omega;
end


% projector for the constraints for A and R
funProjR = @(X) projSimplexR(X,k,p);
funProjA = @(X) projSimplexA(X,qA,EQ_A,n,k,simFlag);


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
fprintf('spectral channels  : %d \n',t);
fprintf('Dictionary size    : %d x %d \n',size(T));
fprintf('epochs             : %d (inner_iterations : %d) \n',epochs,Opt.maxIter);
fprintf('acceleration param : %2.2e \n',rho);
fprintf('Ortho-penalty param: %2.2e \n',lambda);
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
    funObjA = @(X) misfitA(X,F,Y+rho*Q,W,k,lambda,normY);
    Ap      = A;
    A       = minConf_SPG(funObjA,A(:),funProjA,Opt);
    A       = reshape(A,n,k);
    
    % update Q
    res  = (Y - W*(A*F));
    Q    = Q + res;
    
%     figure(2);subplot(2,2,1);plot(F');
%     legend('Orientation','horizontal','Location','bestoutside');
%     subplot(2,2,2);plot(sum(R,2));ylim([0 2]);set(gca,'ytick',[1]);
%     for j=1:k
%         Aj = reshape(A(:,j),[64 64]);
%         subplot(2,k,k+j);imagesc(Aj);axis image;colormap gray;
%         axis off;title(sprintf('%.2f',sum(Aj(:))/sum(A(:))));
%     end
%     pause(0.001);
    
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
    if 1 %(rem(i,floor(epochs/1))==0)
        if printFig, plotFA(F,A,hist,R); end
        fprintf('%10d %12.3e %12.3e %12.3e [%4.2e/%4.2e]\n',i,...
            hist.res(i),hist.res_prog(i),hist.prog(i),time_iter,time_hist);
    end
    
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
        if printFig, plotFA(F,A,hist); end
        fprintf('%10d %12.3e %12.3e %12.3e\n',i,hist.res(i),...
            hist.res_prog(i),hist.prog(i));
        break;
    end
    
    %rho, hist.res(i)
    
    if (rho > 5e-2) && (hist.res_prog(i) > 1e-3)
        rho = 1e-2;
        fprintf('Restoring the original rho value \n');
    end
    if (hist.res_prog(i) < 1e-5)
        rho = rho*10;
        fprintf('Increasing rho by factor 10\n');
    end
    
end

if yWeight, Fsol = Fsol*inv(Omega);end
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

function [f,g] = misfitA(A,F,Y,W,k,lambda,nY)
% misfit function handle for matrix A

% get the size of U and reshape A
[m,n] = size(W);
A     = reshape(A,n,k);

res = W*A*F - Y;                    % residual
f   = 0.5 * norm(res,'fro')^2/nY;   % function value (least-squares misfit)
g   = W' *(res * F');               % gradient wrt A
g   = g(:)/nY;

if lambda > 0
    f1  = sum(sum(A.*(1-A)));
    g1  = 1 - 2*A;

    f   = f + lambda*f1;
    g   = g + lambda*g1(:);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x] = projSimplexA(x,q,EQ,m,n,simplexFlag)
% projector for matrix A (simplex constraints - on rows)

if simplexFlag
    A = reshape(x,m,n);         % reshape matrix R
    B = proj_simplex_q(A',q,EQ);% project onto simplex (transpose required!)
    x = vec(B');                % transpose and vectorize it
else
    x(x < 0) = 0;
    x(x > q) = q;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x] = projSimplexR(x,m,n)

R = reshape(x,m,n);         % reshape matrix R

R = projectTwoSimplexConstr(R);
x = R(:);

% em = ones(m,1);
% en = ones(n,1);
% 
% R  = double(R);
% cvx_begin quiet
%     variable Z(m*n) nonnegative
%     minimize (norm(Z - R(:)))
%     subject to
%         reshape(Z,m,n)*en <= em;
%         reshape(Z,m,n)'*em <= en;
% cvx_end
% x = Z(:);

end

function [x] = projectTwoSimplexConstr(R)

[m,n] = size(R);

MAX_ITER = 1e4;

%%% solve using ADMM

x = R;
if sum(x(:)) >= 1
    for i=1:MAX_ITER
        
        xp=x;
        x1 = transpose(proj_simplex_q(x',1,1));
        x = proj_simplex_q(x1,1,0);
        
        if norm(x1-x,'fro') < 1e-6
            break;
        end
    end
end
%    % update x
%    xp= x;
%    x = proj_simplex_q(transpose((R + rho*(z-u))/(1+rho)),1,0);
%    x = transpose(x);
%    
%    % update z
%    zp= z;
%    z = proj_simplex_q(x+u,1,0);
%    
%    % update u
%    u = u + (x - z);
%    
%    optTol = norm(x - z,'fro');
%    progTol= norm(x-xp,'fro') + norm(z-zp,'fro');
%    
%    if optTol < epTol || progTol < epTol
%        break;
%    end
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


 function p = project_monotone(x)
%function p = project_monotone(x,dir) computes along first direction
%
% This procedure computes the projection onto the constraint set:
%
%                  x(1) <= x(2) <= ... <= x(N)
%
%  INPUTS
% ========
%  x   - ND array

S = cumsum_matrix(x);
S = cummin(S,2,'reverse');
p = max(S,[],1);
p = reshape(p, size(x));

 end

function S = cumsum_matrix(x)
% This function computes the sum:
%
%    S(j,k,:) = sum( x(j:k,:) ) / (k-j+1)
%
% for every 'j' and 'k' in {1,...,N}, where N = size(x,1).


sz = size(x);
S = -Inf([sz(1) sz]);

s = cumsum([zeros([1 sz(2:end)]); x]);
s = reshape(s, [1 sz(1)+1 sz(2:end)]);

for i=1:sz(1)
    S(i,i:end,:) = bsxfun( @minus, s(:,i+1:end,:), s(:,i,:) );
    S(i,i:end,:) = bsxfun( @rdivide, S(i,i:end,:), 1:sz(1)-i+1 );
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plotFA(F,A,hist)

% get sizes
m = size(A,1);
n = size(F,2);
k = size(A,2);

ncol = 6;
nrow = 2*ceil(k/ncol)+1;
Fmax = max(F(:));

% plot
figure(99);m0 = sqrt(m); colormap redblue;
set(gcf,'Position',[696 409 1108 681]);
for i=1:k
    % plot A
    subplot(nrow,ncol,i); Ai = A(:,i); 
    imagesc(reshape(Ai,m0,m0),[-1 1]); 
    axis equal tight;set(gca,'xtick',[],'ytick',[]);
    xlabel(sprintf('A-%d',i));
    
    % plot F
    subplot(nrow,ncol,ncol*ceil(k/ncol)+i);plot(F(i,:),'LineWidth',2); 
    axis([1 n 0 Fmax]);set(gca,'xtick',[],'ytick',[]);
    xlabel(sprintf('F-%d',i));pbaspect([1 1 1])
end
subplot(nrow,ncol,ncol*(nrow-1)+1);semilogy(hist.res,'LineWidth',2);axis square;
title('residual');
subplot(nrow,ncol,ncol*(nrow-1)+2);semilogy(hist.prog,'LineWidth',2);axis square;
title('progress');
subplot(nrow,ncol,ncol*(nrow-1)+3);semilogy(hist.qnorm,'LineWidth',2);axis square;
title('Q');
subplot(nrow,ncol,ncol*(nrow-1)+4);semilogy(hist.anorm,'LineWidth',2);axis square;
title('A');

subplot(nrow,ncol,ncol*(nrow-1)+5);AF1 = A*F(:,1); 
imagesc(reshape(AF1,m0,m0),[-max(AF1) max(AF1)]); 
axis equal tight;set(gca,'xtick',[],'ytick',[]);
% subplot(nrow,ncol,ncol*(nrow-1)+4);semilogy(hist.rnorm,'LineWidth',2);axis square;
% title(sprintf('R-norm (%.4f)',hist.rnorm(end)));
% subplot(nrow,ncol,ncol*(nrow-1)+5);semilogy(hist.rho,'LineWidth',2);axis square;
% title(sprintf('rho (%.4f)',hist.rho(end)));
pause(0.001);

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

% function [x] = projSimplexA2(x,q,EQ,m,n,simplexFlag)
% 
% em= ones(m,1);
% en= ones(n,1);
% 
% cvx_begin quiet
%     variable Z(m*n) nonnegative
%     minimize (norm(Z - x))
%     subject to
%         reshape(Z,m,n)*en == em;
%         reshape(Z,m,n)'*em >= 1;  % floor(0.01*m)*en;
% cvx_end
% 
% x = Z(:);
% 
% end
% 
% function [x] = projSimplexR(x,q,EQ,m,n)
% % project x onto the simplex constraints
% 
% R = reshape(x,m,n);         % reshape matrix R
% B = proj_simplex_q(R',q,EQ);% project onto simplex (transpose required!)
% R = B';                     % transpose it (to reverse back)
% x = R(:);                   % vectorize it
% 
% end

function [X] = thresholdA(X)

X(X < 0.5) = 0;
X(X > 0.5) = 1;

end
