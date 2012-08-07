function [x,itnpqn,f,g,r,xNorm1part,rNorm2part,timePqnProject,testUpdateTau,nLineTot] = minConF_PQN_pqnl1(funObj,x,funProj,options,itn,maxitn,weights,decTol,sigma,b,tau,bpTol,flag)
% 
% Function for using a limited-memory projected quasi-Newton to solve problems of the form
%   min funObj(x) s.t. x in C
%
% The projected quasi-Newton sub-problems are solved the spectral projected
% gradient algorithm
%
%   @funObj(x): function to minimize (returns gradient as second argument)
%   @funProj(x): function that returns projection of x onto C
%
%   options:
%       verbose: level of verbosity (0: no output, 1: final, 2: iter (default), 3:
%       debug)
%       optTol: tolerance used to check for progress (default: 1e-6)
%       maxIter: maximum number of calls to funObj (default: 500)
%       maxProject: maximum number of calls to funProj (default: 100000)
%       numDiff: compute derivatives numerically (0: use user-supplied
%       derivatives (default), 1: use finite differences, 2: use complex
%       differentials)
%       suffDec: sufficient decrease parameter in Armijo condition (default: 1e-4)
%       corrections: number of lbfgs corrections to store (default: 10)
%       adjustStep: use quadratic initialization of line search (default: 0)
%       bbInit: initialize sub-problem with Barzilai-Borwein step (default: 1)
%       SPGoptTol: optimality tolerance for SPG direction finding (default: 1e-6)
%       SPGiters: maximum number of iterations for SPG direction finding (default:10)

nVars = length(x);
timePqnProject = 0;
xNorm1part = zeros(maxitn,1);% why using maxitn instead of opt.maxIter
rNorm2part = zeros(maxitn,1);
itnpqn = 0;% must be initialized here, in case pqn run less than one iter
testUpdateTau = 0;
nLineTot = 0;


% Set Parameters
if nargin < 4
    options = [];
end
[verbose,numDiff,optTol,progTol,maxIter,maxProject,suffDec,corrections,adjustStep,bbInit,SPGoptTol,SPGprogTol,SPGiters,SPGtestOpt] = ...
    myProcessOptions(...
    options,'verbose',2,'numDiff',0,'optTol',1e-6,'progTol',1e-12,'maxIter',100,'maxProject',100000,'suffDec',1e-4,...
    'corrections',inf,'adjustStep',1,'bbInit',1,'SPGoptTol',1e-6,'SPGprogTol',1e-12,'SPGiters',500,'SPGtestOpt',1);

% Output Parameter Settings
if verbose >= 3
   fprintf('Running PQN...\n');
   fprintf('Number of L-BFGS Corrections to store: %d\n',corrections);
   fprintf('Spectral initialization of SPG: %d\n',bbInit);
   fprintf('Maximum number of SPG iterations: %d\n',SPGiters);
   fprintf('SPG optimality tolerance: %.2e\n',SPGoptTol);
   fprintf('SPG progress tolerance: %.2e\n',SPGprogTol);
   fprintf('PQN optimality tolerance: %.2e\n',optTol);
   fprintf('PQN progress tolerance: %.2e\n',progTol);
   fprintf('Quadratic initialization of line search: %d\n',adjustStep);
   fprintf('Maximum number of function evaluations: %d\n',maxIter);
   fprintf('Maximum number of projections: %d\n',maxProject);
end
% Output Log
if verbose >= 2
        fprintf('\n')
        fprintf('Inside of minConf_PQN \n')
        fprintf('%10s %10s %10s %15s %15s %15s\n','Iteration','FunEvals','Projections','Step Length','rNorm2','Opt Cond');
end

% Make objective function (if using numerical derivatives) projects,t,rNorm,optCond
funEvalMultiplier = 1;
if numDiff
    if numDiff == 2
        useComplex = 1;
    else
        useComplex = 0;
    end
    funObj = @(x)autoGrad(x,useComplex,funObj);
    funEvalMultiplier = nVars+1-useComplex;
end

% Project initial parameter vector
x = funProj(x);
projects = 1;

% Evaluate initial parameters
[f,g,r] = funObj(x);
funEvals = 1;

% Check Optimality of Initial Point
projects = projects+1;
if sum(abs(funProj(x-g)-x)) < optTol
    if verbose >= 1
        fprintf('First-Order Optimality Conditions Below optTol at Initial Point\n');
    end
    testUpdateTau = 1;
    return;
end

i = itn;
itnpqn = 0;  % initialize the iteration with 0

while funEvals <= maxIter

        if itnpqn == 0 
            p = funProj(x-g); % dx as spgl1 
            projects = projects+1;
            S = zeros(nVars,0);
            Y = zeros(nVars,0);
            Hdiag = 1;
        else 
            y = g-g_old;
            s = x-x_old;
            [S,Y,Hdiag] = lbfgsUpdate(y,s,corrections,verbose==3,S,Y,Hdiag);
            % Make Compact Representation
            k = size(Y,2);
            L = zeros(k);
            for j = 1:k
                L(j+1:k,j) = S(:,j+1:k)'*Y(:,j);
            end
            N = [S/Hdiag Y];
            M = [S'*S/Hdiag L;L' -diag(diag(S'*Y))];
            HvFunc = @(v)lbfgsHvFunc2(v,Hdiag,N,M);

            if bbInit
                % Use Barzilai-Borwein step to initialize sub-problem
                alpha = (s'*s)/(s'*y);
                if alpha <= 1e-10 || alpha > 1e10
                    alpha = 1/norm(g);
                end

                % Solve Sub-problem
                xSubInit = x-alpha*g;
                feasibleInit = 0;
            else
                xSubInit = x;
                feasibleInit = 1;
            end
            % Solve Sub-problem
            tStart = toc;
            %fprintf('begin solving sub-problem for i = %3.0f\n',i)
            [p,subProjects] = solveSubProblem(x,g,HvFunc,funProj,SPGoptTol,SPGiters,SPGtestOpt,feasibleInit,x,r,decTol,flag);
            projects = projects+subProjects;
            timePqnProject = timePqnProject + toc-tStart;
            %fprintf('finish solving sub-problem for i = %3.0f\n',i)
        end
       
   
        
    d = p-x;
    g_old = g;
    x_old = x;
    
    % Check that Progress can be made along the direction
    gtd = g'*d;
    if gtd >  -progTol
        if verbose >= 1
            fprintf('Directional Derivative below optTol\n');
        end
        %%%%%%%%%orignally I have this here%%%%%%%%%%%%%
%         testUpdateTau = 1;
%         i = i+1;
%         break;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%also, we can go back to previous hession until reach a proper direction 
        if ~isempty(sigma)
            testUpdateTau = 1; 
            i = i+1;
            break;
        else
            while g'*d > -optTol && size(S,2) >= 1
                S = S(:,1:end-1);
                Y = Y(:,1:end-1);
                k = size(Y,2);
                L = zeros(k);
                for j = 1:k
                    L(j+1:k,j) = S(:,j+1:k)'*Y(:,j);
                end
                N = [S/Hdiag Y];
                M = [S'*S/Hdiag L;L' -diag(diag(S'*Y))];
                HvFunc = @(v)lbfgsHvFunc2(v,Hdiag,N,M);
                
                if bbInit
                    % Use Barzilai-Borwein step to initialize sub-problem
                    alpha = (s'*s)/(s'*y);
                    if alpha <= 1e-10 || alpha > 1e10
                        alpha = 1/norm(g);
                    end
                    
                    % Solve Sub-problem
                    xSubInit = x-alpha*g;
                    feasibleInit = 0;
                else
                    xSubInit = x;
                    feasibleInit = 1;
                end
                % Solve Sub-problem
                tStart = toc;
                %fprintf('begin solving sub-problem for i = %3.0f\n',i)
                [p,subProjects] = solveSubProblem(x,g,HvFunc,funProj,SPGoptTol,SPGiters,SPGtestOpt,feasibleInit,x,r,flag);
                projects = projects+subProjects;
                timePqnProject = timePqnProject + toc-tStart;
                
                d = p - x;
            end
            
            % if does not work no matter how far back we go
            if size(S,2) == 1
                d = -g;
                funEvals <= maxIter;
            end

        end




%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       or, quit this round of hessian approxmation
      %d = -g;
      %funEvals = maxIter;
      
      
      
    end
 
    % pqnLineCurvy
   
    f_old = f;
    [f_new,x_new,g_new,r_new,nLine,t,lnErr] = ...
        spgLineCurvy(x,-d,f,funObj,funProj);
    nLineTot = nLineTot + nLine;
    
    % lnErr = 1;
    % curvy line search failed original pqn line search
    if lnErr
             % Select Initial Guess to step length
    if i == itn || adjustStep == 0
       t = 1; 
    else
        t = min(1,2*(f-f_old)/gtd);
    end
    
    % Bound Step length on first iteration
    if i == itn
        t = min(1,1e3/sum(abs(g)));
    end

    % Evaluate the Objective and Gradient at the Initial Step
    x_new = x + t*d;
    [f_new,g_new,r_new] = funObj(x_new);
    funEvals = funEvals+1;

   % Backtracking Line Search
    f_old = f;
     while f_new > f + suffDec*t*gtd || ~isLegal(f_new)
         temp = t;
         nLineTot = nLineTot +1;
         
         % Backtrack to next trial value
         if ~isLegal(f_new) || ~isLegal(g_new)
             if verbose == 3
                 fprintf('Halving Step Size\n');
             end
             t = t/2;
         else
             if verbose == 3
                 fprintf('Cubic Backtracking\n');
             end
             t = polyinterp([0 f gtd; t f_new g_new'*d]);
         end
         
         % Adjust if change is too small/large
         if t < temp*1e-3
             if verbose == 3
                 fprintf('Interpolated value too small, Adjusting\n');
             end
             t = temp*1e-3;
         elseif t > temp*0.6
             if verbose == 3
                 fprintf('Interpolated value too large, Adjusting\n');
             end
             t = temp*0.6;
         end
         
         
         % Check whether step has become too small
         if sum(abs(t*d)) < optTol || t == 0
             if verbose == 3
                 fprintf('Line Search failed\n');
             end
             t = 0;
             f_new = f;
             g_new = g;
             r_new = r;
             break;
         end
         
         % Evaluate New Point
         x_new = x + t*d;
         [f_new,g_new,r_new] = funObj(x_new);
         funEvals = funEvals+1;
         
     end
     end      
    
    % Take Step
    x = x_new;
    f = f_new;
    g = g_new;
    r = r_new;

    dNorm   = norm(d,inf); 
    rNorm   = norm(r, 2);
    bNorm   = norm(b,2);
    gap     = dot(r,(r-b)) + tau*dNorm;
    rGap    = abs(gap) / max(1,f);
    aError1 = rNorm - sigma;
    aError2 = f - sigma^2/2; %%% Lina :: the only difference
    rError1 = abs(aError1) / max(1,rNorm);
    rError2 = abs(aError2) / max(1,f);
    
    optCond = sum(abs(funProj(x-g)-x));
    projects = projects+1;

    i = i + 1;
    itnpqn = itnpqn+1;
    
    % update history
    if isempty(x)
        xNorm1part(i-itn) = 0;
    else
        xNorm1part(i-itn) = NormL1_primal(x,weights);
    end
    rNorm2part(i-itn) = norm(r,2);
    
     % Output Log
    if verbose >= 2
            fprintf('%10d %10d %10d %15.5e %15.5e %15.5e\n',i,funEvals*funEvalMultiplier,projects,t,rNorm,optCond);
    end
    
    % Check optimality
    
    if optCond < optTol
        fprintf('First-Order Optimality Conditions Below optTol\n');
        break;
    end
    
    if sum(abs(t*d)) < progTol
        if verbose >= 1
            fprintf('Step size below optTol\n');
        end
        break;
    end
    
    
    if ~isempty(sigma)
        
        stat = 0;

        if rGap <= max(optTol, rError2) || rError1 <= optTol
            % The problem is nearly optimal for the current tau.
            % Check optimality of the current root.
            test1 = rNorm       <=   bpTol * bNorm;
            test2 = dNorm       <=   bpTol * rNorm;
            test3 = rError1     <=  optTol;
            test4 = rNorm       <=  sigma;

            if test4, stat=1; end % Found suboptimal BP sol.
            if test3, stat=1; end % Found approx root.
            if test2, stat=1; end % Gradient zero -> BP sol.
            if test1, stat=1; end % Resid minim'zd -> BP sol.
        end

        testRelChange1 = (abs(f - f_old) <= decTol * f);
        testRelChange2 = (abs(f - f_old) <= 1e-1 * f * (abs(norm(r,2) - sigma)));
        %testRelChange2 = (abs(f - f_old) <= 1e-1 * f * (abs(rNorm2part(i-itn) - sigma)));
        testUpdateTau  = ((testRelChange1 && rNorm >  2 * sigma) || ...
            (testRelChange2 && rNorm <= 2 * sigma)) && ...
            ~stat;% && ~testUpdateTau;


        if testUpdateTau
            fprintf('% break of testUpdateTau');
            break; % break of testUpdateTau
        end
    end

    if abs(f-f_old) < progTol
        if verbose >= 1
            fprintf('Function value changing by less than optTol\n');
        end
        break;
    end

    if funEvals*funEvalMultiplier > maxIter
        if verbose >= 1
            fprintf('Function Evaluations exceeds maxIter\n');
        end
        break;
    end
    
    if projects > maxProject
        if verbose >= 1
            fprintf('Number of projections exceeds maxProject\n');
        end
        break;
    end
    
    if i >= maxitn;
        break; % break of maximum itns
    end
    
    if norm(r,2) <= optTol
        if verbose >= 1
            fprintf('Optimal solution found\n');
        end
        break;
    end
    
end
% clean up unneeded xnorm1 and rnorm2
itnpqn = i - itn;
if itnpqn == 0 1;
else
    xNorm1part(i-itn+1:end) = [];
    rNorm2part(i-itn+1:end) = [];
end

end



function [p,subProjects] = solveSubProblem(x,g,H,funProj,optTol,maxIter,testOpt,feasibleInit,x_init,r,decTol,flag)
% Uses SPG to solve for projected quasi-Newton direction

options.verbose = 0; 

% if flag
%     options.optTol = max(optTol,1e-3*norm(r)^2);
% else
%     options.optTol = optTol;
% end

options.optTol = optTol;

options.maxIter = maxIter*100;
options.testOpt = testOpt;
options.feasibleInit = feasibleInit;

funObj = @(p)subHv(p,x,g,H);
[p,~,~,subProjects] = minConF_SPG(funObj,x_init,funProj,options,decTol,flag);
end

function [f,g] = subHv(p,x,g,HvFunc)
d = p-x;
Hd = HvFunc(d);
f = g'*d + (1/2)*d'*Hd;
g = g + Hd;
end

function [fNew,xNew,gNew,rNew,iter,step,err] = ...
    spgLineCurvy(x,g,fMax,funObj,funProj,params)
% Projected backtracking linesearch.
% On entry,
% g  is the (possibly scaled) steepest descent direction.

EXIT_CONVERGED  = 0;
EXIT_ITERATIONS = 1;
EXIT_NODESCENT  = 2;
gamma  = 1e-4;
maxIts = 10;
step   =  1;
sNorm  =  0;
scale  =  1;      % Safeguard scaling.  (See below.)
nSafe  =  0;      % No. of safeguarding steps.
iter   =  0;
debug  =  false;  % Set to true to enable log.
n      =  length(x);

if debug
   fprintf(' %5s  %13s  %13s  %13s  %8s\n',...
           'LSits','fNew','step','gts','scale');  
end
   
while 1

    % Evaluate trial point and function value.
    xNew     = funProj(x - step*scale*g);
    [fNew,gNew,rNew] = funObj(xNew);
    s        = xNew - x;
    gts      = scale * real(g' * s);
%   gts      = scale * (g' * s);
    if gts >= 0
       err = EXIT_NODESCENT;
       break
    end

    if debug
       fprintf(' LS %2i  %13.7e  %13.7e  %13.6e  %8.1e\n',...
               iter,fNew,step,gts,scale);
    end
    
    % 03 Aug 07: If gts is complex, then should be looking at -abs(gts).
    % 13 Jul 11: It's enough to use real part of g's (see above).
    if fNew < fMax + gamma*step*gts
%   if fNew < fMax - gamma*step*abs(gts)  % Sufficient descent condition.
       err = EXIT_CONVERGED;
       break
    elseif iter >= maxIts                 % Too many linesearch iterations.
       err = EXIT_ITERATIONS;
       break
    end
    
    % New linesearch iteration.
    iter = iter + 1;
    step = step / 2;

    % Safeguard: If stepMax is huge, then even damped search
    % directions can give exactly the same point after projection.  If
    % we observe this in adjacent iterations, we drastically damp the
    % next search direction.
    % 31 May 07: Damp consecutive safeguarding steps.
    sNormOld  = sNorm;
   sNorm     = norm(s) / sqrt(n);
    %   if sNorm >= sNormOld
    if abs(sNorm - sNormOld) <= 1e-6 * sNorm
       gNorm = norm(g) / sqrt(n);
       scale = sNorm / gNorm / (2^nSafe);
       nSafe = nSafe + 1;
    end
    
end % while 1

end % function spgLineCurvy
