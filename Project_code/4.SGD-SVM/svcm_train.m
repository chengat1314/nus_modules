function [a, b, g, inds, inde, indw] = svcm_train(x, y, C);
% function [a, b, g, inds, inde, indw] = svcm_train(x, y, C);
%        support vector classification machine
%        incremental learning, and leave-one-out cross-validation
%        soft margin
%        uses "kernel.m"
%
%        x: independent variable, (L,N) with L: number of points; N: dimension
%        y: dependent variable, (L,1) containing class labels (-1 or +1)
%        C: soft-margin regularization constant
%
%        a: alpha coefficients (to be multiplied by y)
%        b: offset coefficient
%        g: derivatives (adding one yields margins for each point)
%        inds: indices of support vectors
%        inde: indices of error vectors
%        indw: indices of wrongly classified leave-one-out vectors

%%%%%%%%%% version 1.11; last revised 06/07/2002; send comments to gert@jhu.edu %%%%%%%%%%

%%% GLOBAL VARIABLES:
global doloo online query terse verbose debug memoryhog visualize eps
if isempty(doloo)
    doloo = 1;                        % perform loo at end of training
end
if isempty(online)
    online = 0;                       % take data in the order it is presented
end
if isempty(query)
    query = 0;                        % choose next point based on margin distribution
end
if isempty(terse)
    terse = 0;                        % print out only final results
end
if isempty(verbose)
    verbose = 0;                      % print out details of intermediate results
end
if isempty(debug)
    debug = 0;                        % use only for debugging; slows it down significantly
end
if isempty(memoryhog)
    memoryhog = 0;                    % use more memory; good only if kernel dominates computation
end
if isempty(visualize)
    visualize = 0;                    % record and plot trajectory of coeffs. a, g, and leave-one-out g
end
if isempty(eps)
    eps = 1e-6;                       % margin "margin"; makes Q strictly positive definite
end

%%% END GLOBAL VARIABLES

[L,N] = size(x);
[Ly,Ny] = size(y);
if Ly~=L
    fprintf('svcm_train error: x and y different number of data points (%g/%g)\n\n', L, Ly);
    return
elseif Ny~=1
    fprintf('svcm_train error: y not a single variable (%g)\n\n', Ny);
    return
elseif any(y~=-1&y~=1)
    fprintf('svcm_train error: y takes values different from {-1,+1}\n\n');
    return
end

eps2 = 2*eps/C;
tol = 1e-6;        % tolerance on derivatives at convergence, and their recursive computation

fprintf('Support vector soft-margin classifier with incremental learning\n')
fprintf('  %g training points\n', L)
fprintf('  %g dimensions\n\n', N)

keepe = debug|memoryhog;        % store Qe for error vectors as "kernel cache"
keepr = debug|memoryhog;        % store Qr for recycled support vectors as "kernel cache"
if verbose
    terse = 0;
    if debug
        fprintf('debugging mode (slower, more memory intensive)\n\n')
    elseif memoryhog
        fprintf('memoryhog active (fewer kernel evaluations, more memory intensive)\n\n')
    end
    fprintf('kernel used:\n')
    help kernel
    fprintf('\n')
end

a = zeros(L,1);                       % coefficients, sparse
b = 0;                                % offset
W = 0;                                % energy function
g = -(1+eps)*ones(L,1);               % derivative of energy function

inds = [];                        % indices of support vectors; none initially
inde = [];                        % indices of error vectors; none initially
indo = (L:-1:1)';                 % indices of other vectors; all initially
indr = [];                        % indices of "recycled" other vectors; for memory "caching"
indl = [];                        % indices of leave-one-out vectors (still to be) considered
indw = [];                        % indices of wrongly classified leave-one-out vectors
ls = length(inds);                % number of support vectors;
le = length(inde);                % number of error vectors;
la = ls+le;                       % both
lo = length(indo);                % number of other vectors;
lr = length(indr);                % number of recycled vectors
lw = length(indw);                % number or wrongly classified leave-one-out vectors
processed = zeros(L,1);           % keeps track of which points have been processed
R = Inf;                          % inverse hessian (a(inds) and b only)
Qs = y';                          % extended hessian; (a(inds) plus b, and all vectors)
Qe = [];                          % same, for inde ("cache" for Qs)
Qr = [];                          % same, for indr ("cache" for Qs)
Qc = [];                          % same, for indc (used for gamma and Qs)
if visualize                      % for visualization
    figure(1)
    hold off
    clf
    axis([-0.1*C, 1.1*C, -1.2, 0.2])
    gctraj = [];
    figure(2)
    hold off
    clf
    figure(3)
    hold off
    clf
end

iter = 0;                         % iteration count
memcount = 0;                     % memory usage
kernelcount = 0;                  % kernel computations, counted one "row" (epoch) at a time
training = 1;                     % first do training recursion ...
leaveoneout = 0;                  % ... then do leave-one-out sequence (with retraining)
indc = 0;                         % candidate vector
indco = 0;                        % leave-one-out vector
indso = 0;                        % a recycled support vector; used as buffer
free = a(indo)>0|g(indo)<0;       % free, candidate support or error vector
left = indo(free);                % candidates left
continued = any(left);
while continued                    % check for remaining free points or leave-one-outs to process
    
    % select candidate indc
    indc_prev = indc;
    if online & indc_prev>0
        if query
            processed(indc_prev) = 1;     % record last point in the history log
        else
            processed(1:indc_prev) = 1;   % record last and all preceding points
        end
    end
    if query
%       [gindc, indc] = max(g(left));     % closest to the margin
        [gindc, indc] = min(g(left));     % greedy; worst margin
        indc = left(indc);
    else
        indc = left(length(left));        % take top of the stack, "last-in, first-out"
    end

    % get Qc, row of hessian corresponding to indc (needed for gamma)
    if keepr & lr>0 & ...         % check for match among recycled vectors
        any(find(indr==indc))
            ir = find(indr==indc);                % found, reuse
            Qc = Qr(ir,:);
            indr = indr([1:ir-1,ir+1:lr]);        % ... remove from indr
            Qr = Qr([1:ir-1,ir+1:lr],:);          % ... and Qr
            lr = lr-1;
    elseif indc==indso            % support vector from previous iteration, leftover in memory
        Qc = Qso;
    elseif ls>0 & ...             % check for match among support vectors
        any(find(inds==indc))
            is = find(inds==indc);                % found, reuse
            Qc = Qs(is+1,:);
    elseif keepe & le>0 & ...     % check for match among stored error vectors
        any(find(inde==indc))
            ie = find(inde==indc);                % found, reuse
            Qc = Qe(ie,:);
    elseif indc~=indc_prev        % not (or no longer) available, compute
        xc = x(indc,:);
        yc = y(indc);
        Qc = (yc*y').*kernel(xc,x);
        Qc(indc) = Qc(indc)+eps2;
        kernelcount = kernelcount+1;
    end

    % prepare to increment/decrement z = a(indc)' or y(indc)*b, subject to constraints.
    % move z up when adding indc ((re-)training), down when removing indc (leave-one-out or g>0)
    upc = ~leaveoneout & (g(indc)<=0);
    polc = 2*upc-1;               % polarity of increment in z
    beta = -R*Qs(:,indc);         % change in [b;a(inds)] per change in a(indc)
    if ls>0
        % move z = a(indc)'
        gamma = Qc'+Qs'*beta;     % change in g(:) per change in z = a(indc)'
        z0 = a(indc);             % initial z value
        zlim = max(0,C*polc);     % constraint on a(indc)
    else % ls==0
        % move z = y(indc)*b and keep a(indc) constant; there is no a(:) free to move in inds!
        gamma = y(indc)*Qs';      % change in g(:) per change in z = y(indc)*b
        z0 = y(indc)*b;           % initial z value
        zlim = polc*Inf;          % no constraint on b
    end
    gammac = gamma(indc);
    if gammac<=-tol
        fprintf('\nsvcm_train error: gamma(indc) = %g <= 0 (Q not positive definite)\n\n', gammac)
    elseif gammac==Inf
        fprintf('\nsvcm_train error: gamma(indc) = %g (Q rank deficient)\n\n', gammac)
        return
    end

    % intrinsic limit: g(indc) = 0, where indc becomes support vector
    if ~leaveoneout               % only consider when training indc, not when removing indc
        zlimc = z0-g(indc)'./gammac;
    else                          % leave-indc-out!
        zlimc = polc*Inf;
    end

    % support vector constraints: 0<=a(inds)<=C
    zlims = Inf*polc;                            % by default, immaterial
    if ls>0
        is = find(inds==indc);
        if any(is)                               % leave-indc-out, remove from inds
            zlims = z0;                          % clamp z; no change to variables
        else
            betaa = beta(2:ls+1);                % beta terms corresponding to a(inds)  (not b)
            void = (betaa==0);                   % void zero betaa values ...
            if any(any(~void))
                warning off % suppress div. by 0
                zlims = z0+(C*(betaa*polc>0)-a(inds))./betaa;
                warning on
                zlims(void) = polc*Inf;          % ... which don't enter the constraints
                [zmins, is] = min(zlims*polc,[],1);
                imin = find(zlims==zmins);
                if length(imin)>1
                    [gmax, imax] = max(abs(betaa(imin)),[],1);
                    is = imin(imax);
                end
                zlims = zmins*polc;              % pick tightest constraint
            end
        end
    end

    % error vector constraints: g(inde)<=0
    zlime = Inf*polc;                            % by default, immaterial
    if le>0
        ie = find(inde==indc);
        if any(ie)                               % leave-indc-out, remove from inde
            zlime = z0;                          % clamp z; no change to variables
        else
            gammae = gamma(inde);
            void = (gammae*polc<0)|(gammae==0);  % void g moving down, or zero gamma...
            if any(any(~void))
                warning off % suppress div. by 0
                zlime = z0-g(inde)./gammae;
                warning on
                zlime(void) = polc*Inf;          % ... which don't enter the constraints
                [zmine, ie] = min(zlime*polc,[],1);
                imin = find(zlime==zmine);
                if length(imin)>1
                    [gmax, imax] = max(abs(gammae(imin)),[],1);
                    ie = imin(imax);
                end
                zlime = zmine*polc;              % pick tightest constraint
            end
        end
    end

    % ordinary vector constraints: g(indo)>=0 (only for those that already are)
    zlimo = Inf*polc;                            % by default, immaterial
    if lo>0
        gammao = gamma(indo);
        void = (indo==indc)|(g(indo)<0)|(gammao*polc>0)|(gammao==0);
                                                 % void c, g negative, g moving up, or zero gamma...
        if online
            void = void|~processed(indo);        % ... or, if online, points not seen previously,...
        end
        if any(any(~void))
            warning off % suppress div. by 0
            zlimo = z0-g(indo)./gammao;
            warning on 
            zlimo(void) = polc*Inf;              % ... which don't enter the constraints
            [zmino, io] = min(zlimo*polc,[],1);
            imin = find(zlimo==zmino);
            if length(imin)>1
                [gmax, imax] = max(abs(gammao(imin)),[],1);
                io = imin(imax);
            end
            zlimo = zmino*polc;                  % pick tightest constraint
        end
    end

    % find constraint-satisfying z
    [z,flag] = min([zlim;zlimc;zlims;zlime;zlimo]*polc);
    z = z*polc;
    if (z-z0)*polc<0
        fprintf('\nsvcm_train error: z-z0 of wrong polarity (Q not positive definite)\n\n')
        return
    end

    if verbose & ~leaveoneout & abs(z-z0)<eps
        fprintf('%g*', flag)            % procrastinating iteration!  no progress made
    end

    % update a, b, g and W from z-z0
    if ls>0                             % z = a(indc)
        a(indc) = z;
        b = b+(z-z0)*beta(1);
        a(inds) = a(inds)+(z-z0)*beta(2:ls+1);
        W = W+(z-z0)*(g(indc)'+0.5*(z-z0)*gammac);                % energy
    else                                % z = y(indc)*b
        b = y(indc)*z;
    end
    g = g+(z-z0)*gamma;                 % update g
    iter = iter+1;
    if visualize & ~leaveoneout
        atraj(1:L,iter)=a;              % record trajectory of a(:) over time
        gtraj(1:L,iter)=g;
        ctraj(iter)=indc;
    end 

    % bookkeeping: move elements across indc, inds, inde and indo, and update R and Qs
    converged = (flag<3);               % done with indc; no other changes in inds/inde
    incl_inds = 0;
    if flag==1                          % a(indc) reaches the limits 0 or C, stop moving
        if upc                                  % a(indc)=C, add to inde
            inde = [inde; indc];
            le = le+1;
            if keepe
                Qe = [Qe; Qc];
            end
            a(indc) = C;                        % should be OK, just to avoid round-off
        else % ~upc                                % a(indc)=0, indc stays in (or moves to) indo
            a(indc) = 0;                        % should be OK, just to avoid round-off
        end
    elseif flag==2                      % add indc to support vectors ...
        incl_inds = 1;
        indb = indc;                            % ... store in buffer indb for now
        Qb = Qc;
    elseif flag==3                      % one of support vectors becomes error or other vector
        indso = inds(is);                       % outgoing inds
        Qso = Qs(is+1,:);                       % could be reused later
        free_indc = (indc==indso);              % leave-indc-out: indc is part of inds
        if beta(is+1)*polc<0 | free_indc        % a(indso)=0 or indso=indc, move to indo
            if ~free_indc
                if keepr
                    indr = [indr; indso];       % also recycle into indr for later use
                    Qr = [Qr; Qs(is+1,:)];
                    lr = lr+1;
                end
                if a(indso)>C/2
                    fprintf('svcm_train error: a(indso)=%g; 0 anticipated\n', a(indso));
                end
                a(indso) = 0;                   % should be OK, just to avoid round-off
            end
            g(indso) = 0;                       % same
            indo = [indo; indso];
            lo = lo+1;
        else % beta(is+1)*polc>0 & ~free_indc   % a(indso)=C, move to inde
            if a(indso)<C/2
                fprintf('svcm_train error: a(indso)=%g; C anticipated\n', a(indso));
            end
            if keepe
                Qe = [Qe; Qs(is+1,:)];          % save to memory cache
            end
            a(indso) = C;                       % should be OK, just to avoid round-off
            g(indso) = 0;                       % same
            inde = [inde; indso];
            le = le+1;
        end
        inds = inds([1:is-1,is+1:ls]);          % remove from inds
        stripped = [1:is,is+2:ls+1];            % also from Qs and R ...
        Qs = Qs(stripped,:);
        ls = ls-1;
        if ls > 0
            if R(is+1,is+1)==0
                fprintf('\nsvcm_train error: divide by zero in R contraction\n')
                R(is+1,is+1)=1e-8;
            end
            R = R(stripped,stripped)-R(stripped,is+1)*R(is+1,stripped)/R(is+1,is+1);
        else % no support vectors left
            R = Inf;
        end
    elseif flag==4                       % one of error vectors becomes support/other vector
        indeo = inde(ie);                        % outgoing inde
        if indc==indeo                           % leave-indc-out
            indo = [indo; indeo];                % add inde(ie) to other vectors
            lo = lo+1;
        else
            incl_inds = 1;                       % add inde(ie) to support vectors ...
            indb = indeo;                        % ... store in buffer indb for now
            if keepe
                Qb = Qe(ie,:);                   % recover from Qe cache
            elseif indb==indso
                Qb = Qso;                        % recover from previous outgoing support vector
            else                                 % not in memory either way--- recompute
                Qb = (y(indb)*y').*kernel(x(indb,:),x);
                Qb(indb) = Qb(indb)+eps2;
                kernelcount = kernelcount+1;
            end
        end
        inde = inde([1:ie-1,ie+1:le]);           % remove from inde
        if keepe
            Qe = Qe([1:ie-1,ie+1:le],:);         % remove from Qe
        end
        le = le-1;
    elseif flag==5                       % one of other vectors becomes support vector
        indoo = indo(io);                        % outgoing indo
        incl_inds = 1;                           % add indo(io) to support vectors ...
        indb = indoo;                            % ... store in buffer indb for now
        if keepr & lr>0 & any(find(indr==indb))  % check for match among recycled vectors
            ir = find(indr==indb);                       % found, reuse
            Qb = Qr(ir,:);
            indr = indr([1:ir-1,ir+1:lr]);               % ... remove from indr
            Qr = Qr([1:ir-1,ir+1:lr],:);                 % ... and Qr
            lr = lr-1;
        elseif indb==indso
            Qb = Qso;                            % recover from previous outgoing support vector
        else                                     % not in memory either way--- recompute
            Qb = (y(indb)*y').*kernel(x(indb,:),x);
            Qb(indb) = Qb(indb)+eps2;
            kernelcount = kernelcount+1;
        end
        indo = indo([1:io-1,io+1:lo]);           % remove from indo
        lo = lo-1;
    end

    if incl_inds                        % move buffer indb into support vectors inds
        inds = [inds; indb];                     % move to inds ...
        ls = ls+1;
        Qs = [Qs; Qb];                           % and also Qs and R ...
        if ls==1                                 % compute R directly
            R = [-Qb(indb), y(indb); y(indb), 0];
        else                                     % compute R recursively
            if flag==2                           % from indc; use beta and gamma
                    pivot = gamma(indb);
            else % flag==4                       % from inde; compute beta and pivot
                beta=-R*Qs(1:ls,indb);
                pivot = [beta',1]*Qs(:,indb);
            end
            if pivot<eps2        % should be eps2 when kernel is singular (e.g., linear)
                fprintf('\nsvcm_train error: pivot = %g < %g in R expansion\n\n', pivot, eps2)
                pivot = eps2;
            end
            R = [R,zeros(ls,1);zeros(1,ls+1)]+[beta;1]*[beta',1]/pivot;
        end
    end

    % minor correction in R to avoid numerical instability in recursion when data is near-singular
    Qss = [[0;y(inds)],Qs(:,inds)];
    R = R+R'-R*Qss*R';

    % indc index adjustments (including leave-one-out)
    if converged & (upc|flag==2)                % indc is now part of inds or inde
        i = find(indo==indc);
        indo = indo([1:i-1,i+1:lo]);  % remove indc from indo
        lo = lo-1;
    elseif keepr
        indr = [indr; indc];                    % recycle again into indr for later use
        Qr = [Qr; Qc];
        lr = lr+1;
    end
    if leaveoneout
        indoc = indo(find(indo~=indc));         % indo other than indc
        free = a(indoc)>0|g(indoc)<0;           % candidate support/error vectors in indoc
        if visualize & ~any(free)
            gctraj = [gctraj,[a(indc);g(indc)]];
        end
        converged = converged & ~any(free);
        if converged                            % leave-indc-out reached a(indc)=0, all others settled
            indl = indl(find(indl~=indc));      % remove indc from indl
            if visualize
                figure(1)
                hold on
                if any(gctraj)
                    plot(gctraj(1,:),gctraj(2,:),'ok',gctraj(1,:),gctraj(2,:),'-k')
                end
                gctraj = [];                     % cleanup for next curve
            end
            if g(indc)<-1                       % if leave-indc-out generates an error ...
                indw = [indw; indc];            % ... store its index
                lw = lw+1;                      % ... and increment leave-one-out error count
            end
        end
    end

    % debugging mode: check for consistency; a, R, g and W
    if debug
        f = find(isnan(a));
        if any(f)
           fprintf('svcm_train error: a(%g) = %g\n', [f, a(f)]')
        end
        f = find(isnan(g));
        if any(f)
           fprintf('svcm_train error: g(%g) = %g\n', [f, g(f)]')
        end
        f = find(a>C);
        if any(f)
           fprintf('svcm_train error: a(%g) = %g > C\n', [f, a(f)]')
        end
        f = find(a<0);
        if any(f)
           fprintf('svcm_train error: a(%g) = %g < 0\n', [f, a(f)]')
        end
        f = inde(find(a(inde)>C|a(inde)<C));
        if any(f)
                fprintf('svcm_train error: a(%g) = %g ~= C\n', [f, a(f)]')
        end
        if abs(y'*a)>tol
                fprintf('svcm_train error: y''*a = %g ~= 0 (tol=%g)\n', y'*a, tol);
        end
        Rdiv = max(max(abs(Qss*R-diag(ones(ls+1,1)))));
        if Rdiv>tol
            if flag==2|flag==4    % support vector added
                fprintf('svcm_train error: divergence %g in R expansion (tol=%g; pivot=%g)\n',...
                         Rdiv, tol, pivot);
            elseif flag==3        % support vector removed
                fprintf('svcm_train error: divergence %g in R contraction (tol=%g)\n', Rdiv, tol);
            else                  % no support vector added or removed; strange...
                fprintf('svcm_train error: divergence %g in R ?  (tol=%g)\n', Rdiv, tol);
            end
        end
        if keepe&keepr            % only check g when Qe and Qr are readily available
            greal = Qs'*[b;a(inds)]-(1+eps);
            if le>0
                greal = greal+Qe'*a(inde);
            end
            if lr>0
                greal = greal+Qr'*a(indr);
            end
            if max(abs(greal-g))>tol
                fprintf('svcm_train error: tolerance (tol=%g) exceeded in computation of g\n', tol);
            end
        end
        f = inds(find(abs(g(inds))>tol));
        if any(f)
            fprintf('svcm_train error: g(%g) = %g ~= 0 (tol=%g)\n', [f, g(f), f*0+tol]');
        end
        f = inde(find(g(inde)>tol));
        if any(f)
            fprintf('svcm_train error: g(%g) = %g > 0 (tol=%g)\n', [f, g(f), f*0+tol]');
        end
        Wreal = 0.5*sum((g-b*y-1-eps).*a);        % energy
        if abs(Wreal-W)>tol*abs(W)
            fprintf('svcm_train error: energy W = %g ~= %g (tol=%g)\n', W, Wreal, tol);
        end
        inda = sort([indo;inds;inde]);
        if any(inda~=(1:L)')
            fprintf('svcm_train error: union [indo;inds;inde] does not equate entire set\n');
        end
        if ls~=length(inds)
            fprintf('svcm_train error: miscount in number of support vectors\n');
        end
        if le~=length(inde)
            fprintf('svcm_train error: miscount in number of error vectors\n');
        end
        if lo~=length(indo)
            fprintf('svcm_train error: miscount in number of other vectors\n');
        end
    end

    memcount = max(memcount,ls+keepe*le+keepr*lr+~(keepe&keepr));        % kernel storage

    if verbose
        fprintf('    c: %g (#%g)', y(indc)>0, indc)
        if leaveoneout
            fprintf(', margin: %6.3g, gamma: %6.3g', g(indc)+1, gammac)
        end
        fprintf(' | s: ')
        fprintf('%g', y(inds)>0)
        if ls<6
            fprintf(' (')
            fprintf('#%g', inds)
            fprintf(')')
        end
        fprintf(' | e: 0:%g, 1:%g\n', sum(y(inde)<0), sum(y(inde)>0))
    end

    if converged                        % indc finished; report
        if ~training & ~leaveoneout     % retraining or tracing back; uninteresting
            % nothing
        elseif ~terse & ~leaveoneout
            fprintf(' [%3g%%] iter. %3g: %3g support vectors; %3g error vectors;  energy %g\n', ...
                100-round(100*length(left)/L), iter, ls, le, W)
        elseif ~terse % & leaveoneout
            fprintf(' [%3g%%] iter. %3g: %3g leave-one-out errors', 100-round(100*length(indl)/L), iter, lw)
            fprintf(' (# %g: margin %g)\n', indc, g(indc)+1)
        elseif ~leaveoneout % & terse
            if ls+le>la
                fprintf('+')
            elseif ls+le<la
                fprintf('-')
            else
                fprintf('=')
            end
        else % leaveoneout & terse
            if g(indc)<-1
                fprintf('x')
            else
                fprintf('o')
            end
        end
    la = ls+le;
    end

    % prepare for next iteration, if any
    indoc = indo(find(indo~=indco));       % indo other than indco (leave-one-out index, if active)
    free = a(indoc)>0|g(indoc)<0;          % candidate support/error vectors in indoc
    if any(free)
        left = indoc(free);                % candidates left, keep (re-)training
        leaveoneout = 0;                   % interrupt leave-one-out, if active
    else % ~any(free)                      % done; finish up and (re-)initiate leave-one-out
        if training                        % first time around (not re-training)
            % print out results of svm training
            if terse
                fprintf('\n\n  %g support vectors; %g error vectors; energy %g\n\n', ...
                                 ls, le, W)
            else
                fprintf('\n%4g epoch kernel evaluations (%3g%% of run-time)\n',...
                        kernelcount, round(kernelcount/(ls+le)*100))
                fprintf(  '%4g epoch vectors in memory  (%3g%% of data)\n\n',...
                        memcount, round(memcount/(N+1)*100))
            end
            if visualize

                % plot a trajectory
                figure(2)  
                h=image(atraj/C*length(gray));
                h2=get(h,'Parent');
                set(h2,'YDir','normal')     % 'image' normally reverts the y axis
                colormap(1-gray)
                xlabel('Iteration')
                ylabel('Coefficients \alpha_{\it{i}}')
                print -deps atraj.eps

                % plot g trajectory
                figure(3)
                gmax = max(max(abs(gtraj)));
                h=image((gtraj/gmax+1)/2*length(gray));
                h2=get(h,'Parent');
                set(h2,'YDir','normal')     % 'image' normally reverts the y axis
                grey = gray;                            
                redblue = min(grey,1-grey);             
                redblue(:,1) = redblue(:,1)+1-grey(:,1);
                redblue(:,2) = 2*redblue(:,2);          
                redblue(:,3) = redblue(:,3)+grey(:,3);  
                colormap(redblue)
                xlabel('Iteration')
                ylabel('Coefficients {\it{g}_{\it{i}}}')
                print -deps gtraj.eps

                save traj atraj gtraj ctraj
            end
            if debug                        % store final result to compare with retraining later
                afinal = a;
                gfinal = g;
                Wfinal = W;
            end

            % initiate leave-one-out, if desired
            if doloo
              leaveoneout = 1;
              indl = [inds;inde];           % support and error vectors (others already correct)
              indw = indl(g(indl)<-1);      % leave-one-out errors so far (g<-1) ...
              lw = length(indw);            % ... their number
              indl = indl(g(indl)>=-1);     % remove errors so far from leave-one-out stack
              indco = indl(length(indl));   % pick first leave-one-out index; top of stack ...
              left = indco;                 % ... and let indc=indco (untrain; upc=0)
            else
              continued = 0;
            end

            training = 0;                   % don't ever visit again!
        elseif ~leaveoneout                 % retrained or traced back
            if debug&~any(find(indo==indco))
                % traced back; compare with previously trained results
                if max(abs(a-afinal))>tol
                        fprintf('svcm_train error: final coeffs. a exceed tolerance (tol=%g)\n', tol);
                end
                if max(abs(g-gfinal))>tol
                        fprintf('svcm_train error: final derivatives g exceed tolerance (tol=%g)\n', tol);
                end
                if abs(W-Wfinal)>tol*abs(W)
                        fprintf('svcm_train error: final energy W = %g ~= %g (tol=%g)\n', Wfinal, W, tol);
                end
            end
            if any(indl)
                indco = indl(length(indl)); % leave-one-out index; top of stack
                left = indco;
                leaveoneout = 1;
            else % ~any(indl) % finished all leave-one-outs; summarize, and done!
                continued = 0;
                ltw = sum(g<-1);            % number of training errors
                if terse
                    fprintf('\n\n  %g leave-one-out errors;', lw)
                    fprintf(' %g training errors\n\n', ltw)
                else % ~terse
                    fprintf('\n%4g training points\n', L)
                    fprintf(  '%4g support/error vectors (%g/%g)\n', la, ls, le)
                    fprintf(  '%4g leave-one-out errors  (%3.1f%%)\n', lw, lw/L*100)
                    fprintf(  '%4g training errors       (%3.1f%%)\n', ltw, ltw/L*100)
                    fprintf('\n%4g epoch kernel evaluations (%3g%% of run-time)\n',...
                        kernelcount, round(kernelcount/(ls+le)*100))
                    fprintf(  '%4g epoch vectors in memory  (%3g%% of data)\n\n',...
                        memcount, round(memcount/(N+1)*100))
                end
                if visualize                % plot g trajectory
                    figure(1)
                    hold off
                    xlabel('\alpha_{\it{c}}')
                    ylabel('\it{g_c}')
                    axis([-0.1*C, 1.1*C, -1.2, 0.2])
                    h=line([0,0,C,C],[-1,0,0,-1]);
                    set(h,'Color',[0 0 0])
                    set(h,'LineStyle',':')
                    set(h,'LineWidth',[0.2])
                    h=line([0,C],[-1,-1]);
                    set(h,'Color',[0 0 0])
                    set(h,'LineStyle','--')
                    set(h,'LineWidth',[1.0])
                    print -deps gctraj.eps
                end
            end
            else % leaveoneout              % leave-one-out procedure
            if flag==1                      % finished leave-one-out, now trace back
                % note: already decremented indl and updated indw/lw above
                leaveoneout = 0;
                indco = 0;
            else                            % not done yet with leave-one-out ...
                left = indco;               % ... continue with indc=indco
            end
        end
    end
end
