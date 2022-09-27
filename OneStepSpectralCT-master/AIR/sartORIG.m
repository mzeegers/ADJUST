function [X info restart] = sart(A,b,K,x0,options)
%SART Simultaneous Algebraic Reconstruction Technique (SART) method
%
%   [X info restart] = sart(A,b,K)
%   [X info restart] = sart(A,b,K,x0)
%   [X info restart] = sart(A,b,K,x0,options)
%
% Implements the SART iteration for the linear system Ax = b:
%
%       x^{k+1} = x^k + lambda_k*T*A^T*M*(b-A*x^k)
%
% with T = V^{-1} and M = W^{-1}, where V is a diagonal matrix with row sums 
% of A in the diagonal, and W is a diagonal matrix with the column sums of
% A in the diagonal.
%
% Input:
%   A          m times n matrix.
%   b          m times 1 vector containing the right-hand side.
%   K          Number of iterations. If K is a scalar, then K is the maximum
%              number of iterations and only the last iterate is saved.
%              If K is a vector, then the largest value in K is the maximum
%               number of iterations and only iterates corresponding to the
%              values in K are saved, together with the last iterate.
%              If K is empty then a stopping criterion must be specified.
%   x0         n times 1 starting vector. Default: x0 = 0.
%   options    Struct with the following fields:
%       lambda      The relaxation parameter. If lambda is a scalar then
%                   the corresponding value is used in each iteration. 
%                   If lambda is a string, then it refers to a method to 
%                   determine lambda in each iteration. For this method the
%                   following strings can be specified:
%                       'line'    : lambda is chosen using line search.
%                       'psi1'    : lambda is chosen using the Psi_1-based 
%                                   relaxation method.
%                       'psi1mod' : lambda is chosen using the modified
%                                   Psi_1-based relaxation method.
%                       'psi2'    : lambda is chosen using the Psi_2-based
%                                   relaxation method.
%                       'psi2mod' : lambda is chosen using the modified 
%                                   Psi_2-based relaxation method.
%       stoprule    Struct containing the following information about the
%                   stopping rule:
%                       type = 'none' : (Default) the only stopping rule
%                                       is the maximum number of iterations.
%                              'DP' : The Discrepancy Principle (DP)
%                              'ME' : The Monotone Error rule (ME)
%                              'NCP': Normalized Cumulative Periodogram.
%                       taudelta = product of tau and delta, only needed
%                                  for DP and ME.
%       nonneg      Logical; if true then nonnegativity in enforced in
%                   each iteration.
%       restart     Struct that can contain the first singular value s1
%                   and the diagonals of the matrices M and T.
% Output:
%   X           Matrix containing the saved iterations
%   info        Information vector with 2 elements.
%               info(1) = 0 : stopped by maximum number of iterations
%                         1 : stopped by NCP-rule
%                         2 : stopped by DP-rule
%                         3 : stopped by ME-rule.
%               info(2) = no. of iterations.
%   restart     Struct containing the largest singular value s1, and the
%               diagonals of the matrices M and T.
%
% See also: landweber, cimmino, cav, drop.

% Maria Saxild-Hansen and Per Chr. Hansen, June 11, 2011, DTU Informatics.

% Reference: A. H. Andersen and A. C. Kak, Simultaneous algebraic
% reconstruction technique (SART): A superior implementation of the ART 
% algorithm, Ultrasonic Imaging, 6 (1984), pp. 81-94.

% Check that at least 3 inputs are given.
if nargin < 3
    error('Too few input arguments')
end

[m n] = size(A);

% Check that the sizes of A and b match.
if size(b,1) ~= m || size(b,2) ~= 1
    error('The size of A and b do not match')
end

if nargin < 4
    % Default value for x0.
    x0 = zeros(n,1);
end

% Check if x0 is empty.
if isempty(x0)
    x0 = zeros(n,1);
elseif size(x0,1) ~= n || size(x0,2) ~= 1
    % Check the size og x0.
    error('The size of x0 does not match the problem')
end

AT = A';
rxk = b - A*x0;

% Default values.
if nargin < 5
    % There must be a maximum number of iterations.
    if isempty(K)
        error('No stopping rule specified')
    else
        Knew = sort(K);
        kmax = Knew(end);
        X = zeros(n,length(K));
    end
    
    % Define the matrices W and Vm.
    Apj = full(sum(abs(A),1));
    Aip = full(sum(abs(A),2));
    W = 1./Aip;
    I = (W == Inf);
    W(I) = 0;
    Vm = 1./Apj';
    I = (Vm == Inf);
    Vm(I) = 0;
    
    % If restart is required as output.
    if nargout == 3
        restart.M = W;
        restart.T = Vm;
        restart.s1 = 1;
    end
    
    % Default value of lambda.
    lambda = 1;
    casel = 1;
    
    % Default stopping rule.
    stoprule = 'NO'; 
    k = 0;
    
    % Default there is no nonnegativity projection.
    nonneg = false;
end

% Check the contents of options if present.
if nargin == 5
    
    % Nonnegativity.
    if isfield(options,'nonneg')
        nonneg = options.nonneg;
    else
        nonneg = false;
    end
    
    % Check if W = M is given as input.
    if isfield(options,'restart') && isfield(options.restart,'M')
        W = options.restart.M;
    else
        Aip = full(sum(abs(A),2));
        W = 1./Aip;
        I = (W == Inf);
        W(I) = 0;
    end
    % Check if Vm = T is given as input.
    if isfield(options,'restart') && isfield(options.restart,'T')
        Vm = options.restart.T;
    else
        Apj = full(sum(abs(A),1));
        Vm = 1./Apj';
        I = (Vm == Inf);
        Vm(I) = 0;
    end
    
    % If restart is required as output.
    if nargout == 3
        restart.M = W;
        restart.T = Vm;
        if isfield(options,'restart') && isfield(options.restart,'s1')
            sigma1tildesquare = options.restart.s1^2;
        else
            sigma1tildesquare = 1;
        end
        restart.s1 = sqrt(sigma1tildesquare);
    end
    
    % Stopping rules.
    if isfield(options,'stoprule') && isfield(options.stoprule,'type')
        stoprule = options.stoprule.type;
        % Check that the stoprule is a string.
        if ischar(stoprule)
            if strncmpi(stoprule,'DP',2)
                % DP stopping rule.
                if isfield(options.stoprule,'taudelta')
                    taudelta = options.stoprule.taudelta;
                else
                    error(['The factor taudelta must be specified when '...
                        'using DP'])
                end
                M12 = sqrt(W);
                M12norm = max(M12);
                k = 0;
                
                rk = M12.*rxk;
                nrk = norm(rk);
                
                % Check that the first iteration should be performed.
                if nrk <= taudelta*M12norm
                    info = [2 k];
                    X = x0;
                    return
                end % end the DP-rule.
                
            elseif strncmpi(stoprule,'ME',2)
                % ME stopping rule.
                if isfield(options.stoprule,'taudelta')
                    taudelta = options.stoprule.taudelta;
                else
                    error(['The factor taudelta must be specified when '...
                        'using ME'])
                end
                M12 = sqrt(W);
                M12norm = max(M12);
                k = 0;
                
                rk = M12.*rxk;
                
                if ~isempty(K)
                    K = K + 1;
                end
                
            elseif strncmpi(stoprule,'NC',2)
                % NCP stopping rule.
                dk = inf; 
                q = floor(m/2);
                c_white = (1:q)'./q;
                k = 0;
                
                % Increase the maximum number of iterations by 1.
                if ~isempty(K)
                    K = [K max(K)+1];
                end
                
            elseif strncmpi(stoprule,'NO',2)
                % No stopping rule.
                if isempty(K) 
                    error('No stopping rule specified')
                end
                k = 0;
                
            else
                % Other stopping rules.
                error('The chosen stopping rule is not valid')
            end % end different stopping rules.
            
        else
            error('The stoprule type must be a string')
        end % end stopping rule is string.
        
        % Determine the maximum number of iterations and initialize the
        % output matrix X.
        if isempty(K)
            kmax = inf;
            X = zeros(n,1);
        else
            Knew = sort(K);
            kmax = Knew(end);
            X = zeros(n,length(K));
        end
    else
        % Determine the maximum number of iterations and initialize the
        % output vector X.
        if isempty(K)
            error('No stopping rule specified')
        else
            Knew = sort(K);
            kmax = Knew(end);
            X = zeros(n,length(K));
            stoprule = 'NO';
            k = 0;
            
        end
    end % end stoprule type specified.
    
    % Determine the relaxation parameter lambda:
    if isfield(options,'lambda')
        lambda = options.lambda;
        
        % if lambda is a scalar.
        if ~ischar(lambda)       
            % Convergence check.
            if lambda <= 0 || lambda >= 2
                warning('MATLAB:UnstableRelaxParam',...
                    ['The lambda value is outside the interval (0,2)'])
            end
            casel = 1;
            
        else
            % Calculates the lambda value according to the chosen method.
            if strncmpi(lambda,'line',4)
                % Method: Line search.
                casel = 2;
                
                % Check if the first singular value is known.
                if ~exist('sigma1tildesquare','var')
                    if isfield(options,'restart') && ...
                            isfield(options.restart,'s1')
                        sigma1tildesquare = options.restart.s1^2;
                    else
                        sigma1tildesquare = 1;
                    end
                    
                    % If restart is required as output.
                    if nargout == 3
                        restart.s1 = sqrt(sigma1tildesquare);
                    end
                end
            
            elseif strncmpi(lambda,'psi1',4)
                % Method ENH psi1.
                casel = 3;
                ite = 0;
                
                % Check if the first singular value is known.
                if ~exist('sigma1tildesquare','var')
                    if isfield(options,'restart') && ...
                            isfield(options.restart,'s1')
                        sigma1tildesquare = options.restart.s1^2;
                    else
                        sigma1tildesquare = 1;
                    end
                    
                    % If restart is required as output.
                    if nargout == 3
                        restart.s1 = sqrt(sigma1tildesquare);
                    end
                end
                
                % Precalculate the roots for the case where the stopping
                % iteration is unknown.
                if isempty(K)
                    % Precalculate z for 1000 iterations.
                    z = calczeta(2:999);
                else
                    % Precalculate the roots for the known number of
                    % iterations.
                    z = calczeta(2:max(K)-1);
                end
                
                % Define the values for lambda according to the psi1
                % strategy modified or not.
                if strncmpi(lambda,'psi1mod',7)

                    nu = 2;
                    lambdak = [sqrt(2); sqrt(2); 
                                nu*2*(1-z)]/sigma1tildesquare;
                else

                    lambdak = [sqrt(2); sqrt(2); 
                                2*(1-z)]/sigma1tildesquare;
                end
                
            elseif strncmpi(lambda,'psi2',4)
                % Method: ENH psi2
                casel = 3;
                ite = 0;
                
                % Check if the first singular value is known.
                if ~exist('sigma1tildesquare','var')
                    if isfield(options,'restart') && ...
                            isfield(options.restart,'s1')
                        sigma1tildesquare = options.restart.s1^2;
                    else
                        sigma1tildesquare = 1;
                    end  
                    
                    % If restart is required as output.
                    if nargout == 3
                        restart.s1 = sqrt(sigma1tildesquare);
                    end
                end
                
                % Precalculate the roots for the case where the stopping
                % iteration is unknown.
                if isempty(K)
                    % Precalculate z for 1000 iterations.
                    kk = 2:999;
                    z = calczeta(kk);
                else
                    % Precalculate the roots for the unknown number of
                    % iterations.
                    kk = 2:max(K)-1;
                    z = calczeta(kk);
                end
                
                % Define the values for lambda according to the psi2
                % strategy modified or not.
                if strncmpi(lambda,'psi2mod',7)

                    nu = 1.5;
                    lambdak = [sqrt(2); sqrt(2); 
                         nu*2*(1-z)./((1-z.^(kk')).^2)]/sigma1tildesquare;
                
                else

                    lambdak = [sqrt(2); sqrt(2);
                           2*(1-z)./((1-z.^(kk')).^2)]/sigma1tildesquare;
                end
            else
                error(['The chosen relaxation strategy ''',...
                    options.lambda,''' is not valid for this method.'])  
            end % end check of the lambda strategies.
        end % end check of the class of lambda.
    else    
        % Define a default constant lambda value.
        casel = 1;
        lambda = 1;
        
    end % end if lambda is a field in options.
    
end % end if nargin includes options.

% Initialize the values.
xk = x0;
stop = 0;
l = 0;
klast = 0;

while ~stop
    % Update the iteration number k.
    k = k + 1;
    % Compute the current iteration.

    if casel == 1
        % sart using constant value of lambda.
        xk1 = xk + lambda*(Vm.*(AT*(W.*rxk)));
    elseif casel == 2
        % SART using Line search
        Mrk = W.*rxk;
        ATMrk = AT*Mrk;        
        ATMrkS = sum(ATMrk.^2.*Vm);
        
        lambdak = (rxk'*Mrk)/ATMrkS;
        xk1 = xk + lambdak*(Vm.*ATMrk);
    elseif casel == 3
        % sart using psi1 or psi2
        ite = ite + 1;
        
        % Check if we need to precalculate another 1000 values of lambda.
        if isempty(K)
            if k ~= 1 && rem(k-1,1000) == 0
                ite = 1;
                kk = k:k+999;
                z = calczeta(kk);
                % Define the next lambdak values.
                if strncmpi(lambda,'psi1',4)
                    if strncmpi(lambda,'psi1mod',7)
                        lambdak = (nu*2*(1-z))/sigma1tildesquare;
                    else
                        lambdak = (2*(1-z))/sigma1tildesquare;
                    end
                else
                    if strncmpi(lambda,'psi2mod',7)
                        lambdak = (nu*2*(1-z)./((1-z.^(kk')).^2))/...
                            sigma1tildesquare;
                    else
                        lambdak = (2*(1-z)./((1-z.^(kk')).^2))/...
                            sigma1tildesquare;
                    end
                end
            end
        end
        xk1 = xk + lambdak(ite)*(Vm.*(AT*(W.*rxk)));
    end
    
    % Nonnegativity projection.
    if nonneg, xk1(xk1<0) = 0; end
    
    % New residual.
    rxk1 = b - A*xk1;
    
    % Stopping rule:
    if strncmpi(stoprule,'DP',2)
        % DP stopping rule.
        rk = M12.*rxk1;
        nrk = norm(rk);
        
        if nrk <= taudelta*M12norm || k >= kmax
            stop = 1;
            if k ~= kmax
                info = [2 k];
            else
                info = [0 k];
            end
        end % end the DP-rule.
        
    elseif strncmpi(stoprule,'ME',2)
        % ME stopping rule.
        rk1 = M12.*rxk1;
        nrk = norm(rk);
        dME = rk'*(rk+rk1);
        
        if dME/nrk <= taudelta*M12norm || k >= kmax
            stop = 1;        
            if k ~= kmax
                info = [3 k-1];
            else
                info = [0 k-1];
            end
        else
            rk = rk1;
        end % end the ME-rule.
        
    elseif strncmpi(stoprule,'NC',2)
        % NCP stopping rule.
        rk = rxk1;
        rkh = fft(rk);
        pk = abs(rkh(1:q+1)).^2;
        c = zeros(q,1);
        for index = 1:q
            c(index) = sum(pk(2:index+1))/sum(pk(2:end));
        end
        
        if dk < norm(c-c_white) || k >= kmax
            stop = 1;
            xk1 = xk;
            
            if k ~= kmax
                info = [1 k-1];
            else
                info = [0 k-1];
            end
        else
            dk = norm(c-c_white);
        end % end NCP-rule.
        
    elseif strncmpi(stoprule,'NO',2)
        % No stopping rule.
        if k >= kmax
            stop = 1;
            info = [0 k];
        end
    end % end stoprule type.
    
    % Residual for next iteration.
    rxk = rxk1;
    
    % If the current iteration is requested saved.
    if (~isempty(K) && k == Knew(l+1)) || stop
        l = l + 1;
        % Saves the current iteration.
        if strncmpi(stoprule,'ME',2)
            X(:,l) = xk;
        elseif strncmpi(stoprule,'NC',2)
            if ~(stop && klast == k-1)
                X(:,l) = xk1;
            else
                l = l - 1;
            end
        else
            X(:,l) = xk1;
        end
        klast = k;
    end
    % Update xk.
    xk = xk1;
end

% Save only the required iterations.
X = X(:,1:l);