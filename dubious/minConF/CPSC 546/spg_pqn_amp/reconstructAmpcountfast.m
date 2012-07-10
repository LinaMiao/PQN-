function [xhat,info] = reconstructAmpcountfast(A, y, T, tol, x, verbose)
% RECONSTRUCTAMP recovers a sparse vector x from few linear measurements y.
%
% xhat = reconstructAmp(A, y, T, tol, x, verbose)
%
%   Arguments:
%       A - measurement matrix
%       y - measurements
%       T - max number of iterations (optional)
%       tol - stopping criteria (optional)
%       x - original vector used to print progress of MSE (optional)
%       verbose - print progress optional
%
% Ulugbek Kamilov, LTHC, EPFL, 2010.

% Set some parameters
if ~ exist('T','var')
    T = 500;
end
if ~ exist('tol','var')
    tol = 0.0001;
end
if ~ exist('x','var')
    x = 0;
end
if ~ exist('verbose','var')
    verbose = 0;
end

nmatvec=0;

% Length of the original signal
N = size(A, 2);

% Length of the measurement vector
n = size(A, 1);

% Initial estimate
xhat = zeros(N, 1);
z = y;

%--------------------------------------------------------------------------
% Prepare figure for plotting
%--------------------------------------------------------------------------
if (verbose)
    fprintf(1, 'reconstructAmp: \n');
    h = figure;
end
%--------------------------------------------------------------------------

% Start estimation
for t = 1:T
    % Pre-threshold value
    gamma = xhat + A'*z;

    % Find n-th largest coefficient of gamma
    threshold = largestElement(abs(gamma), n);

    % Estimate the signal (by soft thresholding)
    xhat = eta(gamma, threshold);

    %----------------------------------------------------------------------
    % Compute error and print
    %----------------------------------------------------------------------
    if(verbose)
        err = mean((x-xhat).^2);
        fprintf(1, '[t=%4.d: MSE = %.4f]\n', t, err);
        figure(h);
        plot(t, err, 'b.');
        hold on;
        xlim([1 T]);
        title('Mean-Squared Error of Recovery');
        xlabel('Iteration');
        ylabel('E');
        drawnow;
    end
    %----------------------------------------------------------------------

    % Update the residual
    r = y - A*xhat;
    z = r + (z/n)*sum(etaprime(gamma, threshold));

    nmatvec   = nmatvec+2;
    xNorm1(t) = norm(xhat,1);
    rNorm2(t) = norm(r,2);
        
    % Stopping criteria
    %if(norm(y - A*xhat)/norm(y) < tol)
%        break;
%    end
    info.tau     = norm(xhat,1);
    info.rNorm   = norm(r,2);
    info.iter    = t;
    info.nProdAAt = nmatvec;
    info.xNorm1   = xNorm1;
    info.rNorm2   = rNorm2;
    if ~isempty(find(isnan(xhat))),
        'NaN'
        break;
    end
end