% The file `lobpts.m` is a combination of the following files from
% [Chebfun](https://github.com/chebfun/chebfun/tree/d020ad0f98c1535ed7671ff2b6e702634469daf1):
% 
% jacpts.m
% lobpts.m
% ultrapts.m
% 
% They are covered by the Chebfun license
% https://github.com/chebfun/chebfun/blob/d020ad0f98c1535ed7671ff2b6e702634469daf1/LICENSE.txt
%
% Copyright (c) 2017, The Chancellor, Masters and Scholars of the University 
% of Oxford, and the Chebfun Developers. All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of the University of Oxford nor the names of its 
%       contributors may be used to endorse or promote products derived from 
%       this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR 
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [x, w, v] = lobpts(n, varargin)
%LOBPTS   Gauss-Legendre-Lobatto quadrature nodes and weights.
%  LOBPTS(N) returns N Legendre-Lobatto points X in [-1,1].
%
%  [X, W] = LOBPTS(N) returns also a row vector W of weights for
%  Gauss-Legendre-Lobatto quadrature.
%
%  [X, W, V] = LOBPTS(N) returns additionally a column vector V of weights in
%  the barycentric formula corresponding to the points X. The weights are
%  scaled so that max(abs(V)) = 1.
%
%  In each case, N should be an integer greater than or equal to 2.
%
% See also CHEBPTS, LEGPTS, JACPTS, LEGPOLY, RADAUPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note:
%   The approach used here is to observe that the Gauss-Lobatto points are
%   precisely the roots of (1-x^2)P'_{n-1}(x), and that the roots of P'_{n-1}(x)
%   are the same as the roots of P^(1,1)_{n-2}(x) [NIST, (18.9.15)], which can
%   be obtained for JACPTS. A similar identity [NIST, (18.9.16)] is used for the
%   computation of the quadrature weights from those of JACPTS, and the missing
%   barycentric weights are determined by enforcing the interpolation of f(x) =
%   x or x^2 at x = 0 in the even or odd case respectively.
%
%    x_j = roots of (1-x^2)P'_{n-1}(x)
%    w_j = { 2/(n*(n-1))                        : x_j = -1, 1
%          { 2/(n*(n-1)) * 1/[P_{n-1}(x_j)]^2   : otherwise
%
%   (Note that the weights for n-2 point Gauss-Jacobi with a = b = 1 satisfy 
%    u_j = C/(1-x_j^2)/[d/dx P^(1,1)_{n-2}(x_j)]^2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: Scaled domains?

%% Trivial cases:
if ( n == 1 )
    error('CHEBFUN:lobpts:notSupported', 'N = 1 is not supported.');
elseif ( n == 2 )
    x = [-1 ; 1];
    w = [1, 1];
    v = [-1 ; 1];
    return
elseif ( n == 3 )
    x = [-1 ; 0 ; 1];
    w = [1, 4, 1]/3;
    v = [-.5 ; 1 ; -.5];
    return
end

%% Call JACPTS():
[x, w, v] = jacpts(n - 2, 1, 1, varargin{:});

%% Nodes:
x = [-1 ; x ; 1];

%% Quadrature weights:
w = [-1, w,  1];
w = w./(1-x.^2).';
w([1 end]) = 2/(n*(n - 1));

%% Barycentric weights:
v = v./(1 - x(2:n-1).^2);
v = v/max(abs(v));
if ( mod(n, 2) )
    v1 = -abs(sum(v.*x(2:end-1).^2)/2);
    sgn = 1;
else
    v1 = -abs(sum(v.*x(2:end-1))/2);
    sgn = -1;
end
v = [v1 ; v ; sgn*v1];

end


function [x, w, v] = jacpts(n, a, b, int, meth)
%JACPTS  Gauss-Jacobi quadrature nodes and weights.
%   X = JACPTS(N, ALPHA, BETA) returns the N roots of the degree N Jacobi
%   polynomial with parameters ALPHA and BETA (which must both be > -1)
%   where the Jacobi weight function is w(x) = (1-x)^ALPHA*(1+x)^BETA.
%
%   [X, W] = JACPTS(N, ALPHA, BETA) returns also a row vector W of weights.
%
%   [X, W, V] = JACPTS(N, ALPHA, BETA) returns additionally a column vector V of
%   weights in the barycentric formula corresponding to the points X.
%
%   JACPTS(N, ALPHA, BETA, INTERVAL, METHOD) or JACPTS(N, ALPHA, BETA, METHOD)
%   allows the user to select which method to use.
%    METHOD = 'REC' uses the recurrence relation for the Jacobi polynomials
%     and their derivatives to perform Newton iteration on the WKB approximation
%     to the roots. Default for N < 100.
%    METHOD = 'ASY' uses the Hale-Townsend fast algorithm based upon asymptotic
%     formulae, which is fast and accurate. Default for N >= 100.
%    METHOD = 'GW' uses the traditional Golub-Welsch eigenvalue method,
%     which is maintained mostly for historical reasons.
%
%   [X, W, V] = JACPTS(N, ALPHA, BETA, [A, B]) scales the nodes and weights for
%       the finite interval [A,B].
%
%   The cases ALPHA = BETA = -.5 and ALPHA = BETA = .5 correspond to
%   Gauss-Chebyshev nodes and quadrature, and are treated specially (as a closed
%   form expression for the nodes and weights is available). ALPHA = BETA = 0
%   calls LEGPTS, which is more efficient. The other cases with ALPHA = BETA call
%   ULTRAPTS, which is also faster.
% 
%   When ALPHA ~= BETA and MAX(ALPHA, BETA) > 5 the results may not be accurate. 
%
% See also CHEBPTS, LEGPTS, LOBPTS, RADAUPTS, HERMPTS, LAGPTS, TRIGPTS, and
% ULTRAPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
% 'REC' by Nick Hale, July 2011
% 'ASY' by Nick Hale & Alex Townsend, May 2012 - see [2].
%
% NOTE: The subroutines DO NOT SCALE the weights (with the exception of GW).
% This is done in the main code to avoid duplication.
%
% References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature
%       rules", Math. Comp. 23:221-230, 1969.
%   [2] N. Hale and A. Townsend, "Fast computation of Gauss-Jacobi 
%       quadrature nodes and weights", SISC, 2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults:
interval = [-1, 1];
method = 'default';
method_set = 0;

if ( a <= -1 || b <= -1 )
    error('CHEBFUN:jacpts:sizeAB', 'Alpha and beta must be greater than -1.')
elseif (a~=b && max(a, b) > 5 )
    warning('CHEBFUN:jacpts:largeAB',...
        'ALPHA~=BETA and MAX(ALPHA, BETA) > 5. Results may not be accurate.')
end


% Check inputs:
if ( nargin > 3 )
    if ( nargin == 5 )
        % Calling sequence = JACPTS(N, INTERVAL, METHOD)
        interval = int;
        method = meth;
        method_set = 1;
    elseif ( nargin == 4 )
        if ( ischar(int) )
            % Calling sequence = JACPTS(N, METHOD)
            method = int;
            method_set = true;
        else
            % Calling sequence = JACPTS(N, INTERVAL)
            interval = int;
        end
    end
    validStrings = {'default', 'GW', 'ASY', 'REC'};
    if ( ~any(strcmpi(method, validStrings)) )
        if ( strcmpi(method, 'GLR') )
            error('CHEBFUN:jacpts:glr', ...
                'The GLR algorithm is no longer supported.');
        end
        error('CHEBFUN:jacpts:inputs', ['Unrecognised input string: ', method]);
    end
end

if ( any(isinf(interval)) )  % Inf intervals not yet supported.
    % TODO: How do we scale the weights?
    error('CHEBFUN:jacpts:infinterval', ... 
          'JACPTS() does not yet support infinite intervals');
elseif ( numel(interval) > 2 )
    warning('CHEBFUN:legpts:domain',...
        'Piecewise intervals are not supported and will be ignored.');
    interval = interval([1, end]);
end

% Deal with trivial cases:
if ( n < 0 )
    error('CHEBFUN:jacpts:n', 'First input should be a positive number.');
elseif ( n == 0 )   % Return empty vectors if n == 0:
    x = []; 
    w = []; 
    v = []; 
    return
elseif ( n == 1 )
    x0 = (b-a)/(a+b+2);
    x = diff(interval)/2 * (x0+1) + interval(1); % map from [-1,1] to interval. 
    w = 2^(a+b+1)*beta(a+1, b+1) * diff(interval)/2;
    v = 1;
    return
end

% Special cases:
if ( a == b )
    if ( a == 0 )  % Gauss-Legendre: alpha = beta = 0
        [x, w, v] = legpts(n, method);
        [x, w] = rescale_jacpts(x, w, interval, a, b);
        return
    elseif ( a == -.5 )  % Gauss-Chebyshev: alpha = beta = -.5
        [x, ignored, v] = chebpts(n, interval, 1);
        w = repmat(pi/n,1,n);
        [ignored, w] = rescale_jacpts(x, w, interval, a, b);
        return
    elseif ( a == .5 )   % Gauss-Chebyshev2: alpha = beta = .5
        x = chebpts(n+2, 2);
        x = x(2:n+1);
        w = pi/(n+1)*(1-x.^2).';
        v = (1-x.^2);
        v(2:2:end) = -v(2:2:end);
        [x, w] = rescale_jacpts(x,w,interval,a,b);
        return
    else % Gauss-Gegenbauer: alpha = beta
        lambda = a + .5;
        [x, w, v] = ultrapts(n, lambda, interval);
        return
    end
end

% Choose an algorithm:
if ( n < 20 || (n < 100 && ~method_set) || strcmpi(method, 'rec_jacpts') )
    [x, w, v] = rec_jacpts(n, a, b); % REC (Recurrence relation)
    
elseif ( strcmpi(method, 'GW') )
    [x, w, v] = gw_jacpts(n, a, b);  % GW  see [1]
    
else
    [x, w, v] = asy_jacpts(n, a, b); % HT  see [2]
    
end

% Compute the constant for the weights:
if ( ~strcmpi(method,'GW') )
    if ( n >= 100 )
        cte1 = ratioOfGammaFunctions_jacpts(n+a+1, b);
        cte2 = ratioOfGammaFunctions_jacpts(n+1, b);
        C = 2^(a+b+1)*(cte2/cte1);
    else
        C = 2^(a+b+1) * exp( gammaln(n+a+1) - gammaln(n+a+b+1) + ...
            gammaln(n+b+1) - gammaln(n+1) );   
        % An alternative approach could be used: 
        %   C = 2^(a+b+1)*beta(a+1, b+1)/sum(w), 
        % but we prefer compute the weights independently.
    end
    w = C*w;
end

% Scale the nodes and quadrature weights:
[x, w] = rescale_jacpts(x, w, interval, a, b);

% Scale the barycentric weights:
v = abs(v); 
v(2:2:end) = -v(2:2:end);
v = v./max(abs(v)); 

end

function [x, w] = rescale_jacpts(x, w, interval, a, b)
%RESCALE   Rescale nodes and weights to an arbitrary finite interval.
    if ( ~all(interval == [-1, 1]) )
        c1 = .5*sum(interval); 
        c2 = .5*diff(interval);
        w = c2^(a+b+1)*w;
        x = c1 + c2*x;    
    end
end

function cte = ratioOfGammaFunctions_jacpts(m,delta)
%RATIOGAMMA Compute the ratio gamma(m+delta)/gamma(m). See [2].
    % cte = gamma(m+delta)/gamma(m)
    ds = .5*delta^2/(m-1);
    s = ds;
    j = 1;
    while ( abs(ds/s) > eps/100 ) % Taylor series in expansion 
        j = j+1;
        ds = -delta*(j-1)/(j+1)/(m-1)*ds;
        s = s + ds;
    end
    p2 = exp(s)*sqrt(1+delta/(m-1))*(m-1)^(delta);
    % Stirling's series:
    g = [1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880, ...
        5246819/75246796800, -534703531/902961561600, ...
        -4483131259/86684309913600, 432261921612371/514904800886784000];
    f = @(z) sum(g.*[1, cumprod(ones(1, 9)./z)]);
    cte = p2*(f(m+delta-1)/f(m-1));
end


%% ------------------------- Routines for GW ----------------------------
    
function [x, w, v] = gw_jacpts(n, a, b)
    ab = a + b;
    ii = (2:n-1)';
    abi = 2*ii + ab;
    aa = [(b - a)/(2 + ab)
          (b^2 - a^2)./((abi - 2).*abi)
          (b^2 - a^2)./((2*n - 2+ab).*(2*n+ab))];
    bb = [2*sqrt( (1 + a)*(1 + b)/(ab + 3))/(ab + 2) ; 
          2*sqrt(ii.*(ii + a).*(ii + b).*(ii + ab)./(abi.^2 - 1))./abi];
    TT = diag(bb,1) + diag(aa) + diag(bb,-1); % Jacobi matrix.
    [V, x] = eig( TT );                       % Eigenvalue decomposition.
    x = diag(x);                              % Jacobi points.
    w = V(1,:).^2*2^(ab+1)*beta(a+1, b+1);    % Quadrature weights.
    v = sqrt(1-x.^2).*abs(V(1,:))';           % Barycentric weights.
end

%% ------------------------- Routines for REC ---------------------------


function [x, w, v] = rec_jacpts(n, a, b)
%REC   Compute nodes and weights using recurrrence relation.

   [x1, ders1] = rec_main_jacpts(n, a, b, 1); % Nodes and P_n'(x)
   [x2, ders2] = rec_main_jacpts(n, b, a, 0); % Nodes and P_n'(x)
   x = [-x2(end:-1:1) ; x1];
   ders = [ders2(end:-1:1) ; ders1];
   w = 1./((1-x.^2).*ders.^2)';        % Quadrature weights
   v = 1./ders;                        % Barycentric weights
   
end

function [x, PP] = rec_main_jacpts(n, a, b, flag)
%REC_MAIN   Jacobi polynomial recurrence relation.

% Asymptotic formula (WKB) - only positive x.
if ( flag )
    r = ceil(n/2):-1:1;
else
    r = floor(n/2):-1:1;  
end
C = (2*r+a-.5)*pi/(2*n+a+b+1);
T = C + 1/(2*n+a+b+1)^2 * ((.25-a^2)*cot(.5*C) - (.25-b^2)*tan(.5*C));
x = cos(T).';

% Initialise:
dx = inf; 
l = 0;
% Loop until convergence:
while ( (norm(dx,inf) > sqrt(eps)/1000) && (l < 10) )
    l = l + 1;
    [P, PP] = eval_Jac_jacpts(x, n, a, b);
    dx = -P./PP; 
    x = x + dx;
end
% Once more for derivatives:
[ignored, PP] = eval_Jac_jacpts(x, n, a, b);

end

function [P, Pp] = eval_Jac_jacpts(x, n, a, b)
%EVALJAC   Evaluate Jacobi polynomial and derivative via recurrence relation.

% Initialise:
ab = a + b;
P = .5*(a-b+(ab+2)*x);  
Pm1 = 1; 
Pp = .5*(ab+2);         
Ppm1 = 0; 

% n = 0 case:
if ( n == 0 )
    P = Pm1; 
    Pp = Ppm1; 
end

for k = 1:n-1
    % Useful values:
    A = 2*(k + 1)*(k + ab + 1)*(2*k + ab);
    B = (2*k + ab + 1)*(a^2 - b^2);
    C = prod(2*k + ab + (0:2)');
    D = 2*(k + a)*(k + b)*(2*k + ab + 2);

    % Recurrence:
    Pa1 = ( (B+C*x).*P - D*Pm1 ) / A;
    Ppa1 = ( (B+C*x).*Pp + C*P - D*Ppm1 ) / A;

    % Update:
    Pm1 = P; 
    P = Pa1;  
    Ppm1 =  Pp; 
    Pp = Ppa1;
end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------- Routines for ASY algorithm ------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v] = asy_jacpts(n, a, b)
%ASY   Compute nodes and weights using asymptotic formulae.

    if ( n <= 20 ) % Use only boundary formula:
        [xbdy, wbdy, vbdy] = asy2_jacpts(n, a, b, ceil(n/2));  
        [xbdy2, wbdy2, vbdy2] = asy2_jacpts(n, b, a, floor(n/2));  
        x = [-xbdy2(end:-1:1) ; xbdy];
        w = [wbdy2(end:-1:1), wbdy];
        v = [vbdy2(end:-1:1) ; vbdy];
        return
    end

    % Determine switch between interior and boundary regions:
    nbdy = 10;
    bdyidx1 = n-(nbdy-1):n; 
    bdyidx2 = nbdy:-1:1;

    % Interior formula:
    [x, w, v] = asy1_jacpts(n,a,b,nbdy);   

    % Boundary formula (right):
    [xbdy, wbdy, vbdy] = asy2_jacpts(n, a, b, nbdy);  
    x(bdyidx1) = xbdy;  
    w(bdyidx1) = wbdy; 
    v(bdyidx1) = vbdy;
    
    % Boundary formula (left):
    if ( a ~= b )
        [xbdy, wbdy, vbdy] = asy2_jacpts(n, b, a, nbdy);  
    end
    x(bdyidx2) = -xbdy; 
    w(bdyidx2) = wbdy; 
    v(bdyidx2) = vbdy;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             ASY1 (Interior)                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v] = asy1_jacpts(n, a, b, nbdy)
% Algorithm for computing nodes and weights in the interior.

    % Approximate roots via asymptotic formula: (Gatteschi and Pittaluga, 1985)
    K = (2*(n:-1:1)+a-.5)*pi/(2*n+a+b+1);
    tt = K + 1/(2*n+a+b+1)^2*((.25-a^2)*cot(.5*K)-(.25-b^2)*tan(.5*K));

    % First half (x > 0):
    t = tt(tt <= pi/2);
    mint = t(end-nbdy+1);
    idx = 1:max(find(t < mint,1)-1, 1);

    dt = inf; j = 0;
    % Newton iteration
    while ( norm(dt,inf) > sqrt(eps)/100 && j < 10 )
        [vals, ders] = feval_asy1_jacpts(n, a, b, t, idx, 0);  % Evaluate
        dt = vals./ders;                                % Newton update
        t = t + dt;                                     % Next iterate
        j = j + 1;
        dt = dt(idx);
    end
    [vals, ders] = feval_asy1_jacpts(n, a, b, t, idx, 1);      % Once more for luck
    t = t + vals./ders;                                 % Newton update.

    % Store:
    x = cos(t);
    w = 1./ders.^2;
    v = (sin(t)./ders);

    % Second half (x < 0):
    tmp = a; 
    a = b; 
    b = tmp;
    t = pi - tt(1:(n-length(x)));
    mint = t(nbdy);
    idx = max(find(t > mint, 1), 1):numel(t);

    dt = inf; j = 0;
    % Newton iteration
    while ( norm(dt,inf) > sqrt(eps)/100 && j < 10 )
        [vals, ders] = feval_asy1_jacpts(n, a, b, t, idx, 0);  % Evaluate.
        dt = vals./ders;                                % Newton update.
        t = t + dt;                                     % Next iterate.
        j = j + 1;
        dt = dt(idx);
    end
    [vals, ders] = feval_asy1_jacpts(n, a, b, t, idx, 1);      % Once more for luck.
    t = t + vals./ders;                                 % Newton update.

    % Store:
    x = [-cos(t) x].';
    w = [1./ders.^2 w];
    v = [sin(t)./ders v].';

end

% -------------------------------------------------------------------------

function [vals, ders] = feval_asy1_jacpts(n, a, b, t, idx, flag)
% Evaluate the interior asymptotic formula at x = cos(t).
    
    % Number of terms in the expansion:
    M = 20;

    % Some often used vectors/matrices:
    onesT = ones(1,length(t));
    onesM = ones(M,1);
    MM = transpose(0:M-1);

    % The sine and cosine terms:
    alpha = (.5*(2*n+a+b+1+MM))*onesT .* (onesM*t) - .5*(a+.5)*pi;
    cosA = cos(alpha);
    sinA = sin(alpha);

    if ( flag ) % Evaluate cos(alpha) using Taylor series.
        k = 1:numel(t);
        if ( idx(1) == 1 )
            k = fliplr(k);
        end
        % Hi-lo computation to squeeze an extra digit in the computation.
        ta = double(single(t));    tb = t - ta;
        hi = n*ta;                 lo = n*tb + (a+b+1)*.5*t; 
        pia = 3.1415927410125732;  pib = -8.742278000372485e-08; % pib = pi-pia;
        % Doing this means that pi - pia - pib = 10^-24
        dh = ( hi - (k-.25)*pia ) + lo - .5*a*pia - (k-.25+.5*a)*pib;
        tmp = 0; sgn = 1; fact = 1; DH = dh; dh2 = dh.*dh;       % Initialise.
        for j = 0:20
            dc = sgn*DH/fact;
            tmp = tmp + dc;
            sgn = -sgn;
            fact = fact*(2*j+3)*(2*j+2);
            DH = DH.*dh2;
            if ( norm(dc, inf) < eps/2000 )
                break
            end
        end
        tmp(2:2:end) = -tmp(2:2:end);          % }
        [~, loc] = max(abs(tmp));              %  } Fix up the sign.
        tmp = sign(cosA(1,loc)*tmp(loc))*tmp;  % }
        cosA(1,:) = tmp;
    end

    sinT = onesM*sin(t);
    cosT = onesM*cos(t);
    cosA2 = cosA.*cosT + sinA.*sinT;
    sinA2 = sinA.*cosT - cosA.*sinT;

    one = ones(1,length(t));
    sinT = [one ; cumprod(onesM(2:end)*(.5*csc(.5*t)))];
    cosT = .5*sec(.5*t);

    j = 0:M-2;
    vec = (.5+a+j).*(.5-a+j)./(j+1)./(2*n+a+b+j+2);
    P1 = [1  cumprod(vec)];
    P1(3:4:end) = -P1(3:4:end);
    P1(4:4:end) = -P1(4:4:end);
    P2 = eye(M);
    for l = 0:M-1
        j = 0:(M-l-2);
        vec = (.5+b+j).*(.5-b+j)./(j+1)./(2*n+a+b+j+l+2);
        P2(l+1+(1:length(j)),l+1) = cumprod(vec);
    end
    PHI = repmat(P1,M,1).*P2;

    j = 0:M-2;
    vec = (.5+a+j).*(.5-a+j)./(j+1)./(2*(n-1)+a+b+j+2);
    P1 = [1  cumprod(vec)];
    P1(3:4:end) = -P1(3:4:end);
    P1(4:4:end) = -P1(4:4:end);
    P2 = eye(M);
    for l = 0:M-1
        j = 0:(M-l-2);
        vec = (.5+b+j).*(.5-b+j)./(j+1)./(2*(n-1)+a+b+j+l+2);
        P2(l+1+(1:length(j)),l+1) = cumprod(vec);
    end
    PHI2 = repmat(P1,M,1).*P2;

    S = 0; S2 = 0;
    SC = sinT;
    for m = 0:M-1

        l = 0:2:m;
        phi = PHI(m+1,l+1);
        dS1 = phi*SC(l+1,:).*cosA(m+1,:);

        phi2 = PHI2(m+1,l+1);
        dS12 = phi2*SC(l+1,:).*cosA2(m+1,:);

        l = 1:2:m;
        phi = PHI(m+1,l+1);
        dS2 = phi*SC(l+1,:).*sinA(m+1,:);

        phi2 = PHI2(m+1,l+1);
        dS22 = phi2*SC(l+1,:).*sinA2(m+1,:);

        if m > 10 && norm(dS1(idx) + dS2(idx),inf) < eps/100, break, end

        S = S + dS1 + dS2;
        S2 = S2 + dS12 + dS22;

        SC(1:m+1,:) = bsxfun(@times,SC(1:m+1,:),cosT);
    end

    % Constant out the front:
    dsa = .5*(a^2)/n; dsb = .5*(b^2)/n; dsab = .25*(a+b)^2/n;
    ds = dsa + dsb - dsab; s = ds; j = 1; 
    dsold = ds; % to fix a = -b bug.
    while ( (abs(ds/s) + dsold) > eps/10 )
        dsold = abs(ds/s);
        j = j+1;
        tmp = -(j-1)/(j+1)/n;
        dsa = tmp*dsa*a;
        dsb = tmp*dsb*b;
        dsab = .5*tmp*dsab*(a+b);
        ds = dsa + dsb - dsab;
        s = s + ds;
    end
    p2 = exp(s)*sqrt(2*pi)*sqrt((n+a)*(n+b)/(2*n+a+b))/(2*n+a+b+1);
    g = [1 1/12 1/288 -139/51840 -571/2488320 163879/209018880 ...
         5246819/75246796800 -534703531/902961561600 ...
         -4483131259/86684309913600 432261921612371/514904800886784000];
    f = @(z) sum(g.*[1 cumprod(ones(1,9)./z)]);
    C = p2*(f(n+a)*f(n+b)/f(2*n+a+b))*2/pi;
    C2 = C*(a+b+2*n).*(a+b+1+2*n)./(4*(a+n).*(b+n));

    vals = C*S;
    S2 = C2*S2;

    % Use relation for derivative:
    ders = (n*(a-b-(2*n+a+b)*cos(t)).*vals + 2*(n+a)*(n+b)*S2)/(2*n+a+b)./sin(t);
    denom = 1./real(sin(t/2).^(a+.5).*cos(t/2).^(b+.5));
    vals = vals.*denom;
    ders = ders.*denom;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             ASY2 (Boundary)                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v] = asy2_jacpts(n, a, b, npts)
% Algorithm for computing nodes and weights near the boundary.

% Use Newton iterations to find the first few Bessel roots:
smallK = min(30, npts);
jk = besselroots(a, min(npts, smallK));
% Use asy_jacpts formula for larger ones (See NIST 10.21.19, Olver 1974 p247)
if ( npts > smallK )
    mu = 4*a^2;
    a8 = 8*((length(jk)+1:npts).'+.5*a-.25)*pi;
    jk2 = .125*a8-(mu-1)./a8 - 4*(mu-1)*(7*mu-31)/3./a8.^3 - ...
          32*(mu-1)*(83*mu.^2-982*mu+3779)/15./a8.^5 - ...
          64*(mu-1)*(6949*mu^3-153855*mu^2+1585743*mu-6277237)/105./a8.^7;
    jk = [jk ; jk2];
end
jk = real(jk(1:npts));

% Approximate roots via asymptotic formula: (see Olver 1974, NIST, 18.16.8)
rho = n + .5*(a + b + 1); 
phik = jk/rho;
t = phik + ((a^2-.25)*(1-phik.*cot(phik))./(2*phik) - ...
    .25*(a^2-b^2)*tan(.5*phik))/rho^2;

% Only first half (x > 0):
if ( any(t > 1.1*pi/2) )
    warning('CHEBFUN:jacpts:asy2_jacpts:theta', ...
        'Theta > pi/2. Result may be inaccurate.'); 
end

% Compute higher terms:
[tB1, A2, tB2, A3] = asy2_higherterms_jacpts(a, b, t, n);

dt = inf; j = 0;
% Newton iteration:
while ( norm(dt,inf) > sqrt(eps)/200 && j < 10)
    [vals, ders] = feval_asy2_jacpts(n, t, 1); % Evaluate via asymptotic formula.
    dt = vals./ders;                    % Newton update.
    t = t + dt;                         % Next iterate.
    j = j + 1;
end
[vals, ders] = feval_asy2_jacpts(n, t, 1);     % Evaluate via asymptotic formula.
dt = vals./ders;                        % Newton update
t = t + dt;    
    
% flip:
t = t(npts:-1:1); 
ders = ders(npts:-1:1);
% vals = vals(npts:-1:1);

% Revert to x-space:
x = cos(t);      
w = (1./ders.^2).';   
v = sin(t)./ders;

    function [vals, ders] = feval_asy2_jacpts(n, t, flag)
    % Evaluate the boundary asymptotic formula at x = cos(t).
    
        % Useful constants:
        rho2 = n + .5*(a + b - 1);
        A = (.25 - a^2);       
        B = (.25 - b^2);
        
        % Evaluate the Bessel functions:
        Ja = besselj(a, rho*t, 0);
        Jb = besselj(a + 1, rho*t, 0);
        Jbb = besselj(a + 1, rho2*t, 0);
        if ( ~flag )
            Jab = besselj(a, rho2*t, 0);
        else
            % In the final step, perform accurate evaluation
            Jab = besselTaylor_jacpts(-t, rho*t, a);
        end

        % Evaluate functions for recurrsive definition of coefficients:
        gt = A*(cot(t/2)-(2./t)) - B*tan(t/2);
        gtdx = A*(2./t.^2-.5*csc(t/2).^2) - .5*B*sec(t/2).^2;
        tB0 = .25*gt;
        A10 = a*(A+3*B)/24;
        A1 = gtdx/8 - (1+2*a)/8*gt./t - gt.^2/32 - A10;
        % Higher terms:
        tB1t = tB1(t); 
        A2t = A2(t); 

        % VALS:
        vals = Ja + Jb.*tB0/rho + Ja.*A1/rho^2 + Jb.*tB1t/rho^3 + Ja.*A2t/rho^4;
        % DERS:
        vals2 = Jab + Jbb.*tB0/rho2 + Jab.*A1/rho2^2 + Jbb.*tB1t/rho2^3 + Jab.*A2t/rho2^4;
        
        % Higher terms (not needed for n > 1000).
        tB2t = tB2(t); A3t = A3(t);
        vals = vals + Jb.*tB2t/rho^5 + Ja.*A3t/rho^6;
        vals2 = vals2 + Jbb.*tB2t/rho2^5 + Jab.*A3t/rho2^6;
        
        % Constant out the front (Computed accurately!)
        ds = .5*(a^2)/n;
        s = ds; jj = 1;
        while abs(ds/s) > eps/10
            jj = jj+1;
            ds = -(jj-1)/(jj+1)/n*(ds*a);
            s = s + ds;
        end
        p2 = exp(s)*sqrt((n+a)/n)*(n/rho)^a;
        g = [1 1/12 1/288 -139/51840 -571/2488320 163879/209018880 ...
             5246819/75246796800 -534703531/902961561600 ...
             -4483131259/86684309913600 432261921612371/514904800886784000];
        f = @(z) sum(g.*[1 cumprod(ones(1,9)./z)]);
        C = p2*(f(n+a)/f(n))/sqrt(2);

        % Scaling:
        valstmp = C*vals;
        denom = sin(t/2).^(a+.5).*cos(t/2).^(b+.5);
        vals = sqrt(t).*valstmp./denom;

        % Relation for derivative:
        C2 = C*n/(n+a)*(rho/rho2)^a;
        ders = (n*(a-b-(2*n+a+b)*cos(t)).*valstmp + 2*(n+a)*(n+b)*C2*vals2)/(2*n+a+b);
        ders = ders.*(sqrt(t)./(denom.*sin(t)));
        
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [TODO]: The codes below are only here temporarily.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ja = besselTaylor_jacpts(t, z, a)
%BESSELTAYLOR    Accurate evaluation of Bessel function J_A for asy2_jacpts. (See [2].)
% BESSELTAYLOR(T, Z, A) evaluates J_A(Z+T) by a Taylor series expansion about Z. 

npts = numel(t);
kmax = min(ceil(abs(log(eps)/log(norm(t, inf)))), 30);
H = bsxfun(@power, t, 0:kmax).';
% Compute coeffs in Taylor expansions about z (See NIST 10.6.7)
[nu, JK] = meshgrid(-kmax:kmax, z);
Bjk = besselj(a + nu, JK, 0);
nck = abs(pascal(floor(1.25*kmax), 1)); nck(1,:) = []; % nchoosek
AA = [Bjk(:,kmax+1), zeros(npts, kmax)];
fact = 1;
for k = 1:kmax
    sgn = 1;
    for l = 0:k
        AA(:,k+1) = AA(:,k+1) + sgn*nck(k,l+1)*Bjk(:,kmax+2*l-k+1);
        sgn = -sgn;
    end
    fact = k*fact;
    AA(:,k+1) = AA(:,k+1)/2^k/fact;
end
% Evaluate Taylor series:
Ja = zeros(npts, 1);
for k = 1:npts
    Ja(k,1) = AA(k,:)*H(:,k);
end
end

function [tB1, A2, tB2, A3, tB3, A4] = asy2_higherterms_jacpts(alph, bet, theta, n)
%ASY2_HIGHERTERMS   Higher-order terms for boundary asymptotic series.
% Compute the higher order terms in asy2_jacpts boundary formula. See [2]. 

% These constants are more useful than alph and bet:
A = (0.25 - alph^2);
B = (0.25 - bet^2);

% For now, just work on half of the domain:
c = max(max(theta), .5);
if ( n < 30 )
    N = ceil(40 - n);
elseif ( n >= 30 && c > pi/2-.5)
    N = 15;
else
    N = 10;
end
Nm1 = N - 1;

% Scaled 2nd-kind Chebyshev points and barycentric weights:
t = .5*c*( sin(pi*(-Nm1:2:Nm1)/(2*Nm1)).' + 1 );
v = [.5 ; ones(Nm1,1)];
v(2:2:end) = -1;
v(end) = .5*v(end);

% The g's:
g = A*(cot(t/2)  -2./t) - B*tan(t/2);
gp = A*(2./t.^2 - .5*csc(t/2).^2) - .5*(.25-bet^2)*sec(t/2).^2;
gpp = A*(-4./t.^3 + .25*sin(t).*csc(t/2).^4) - 4*B*sin(t/2).^4.*csc(t).^3;
g(1) = 0; gp(1) = -A/6-.5*B; gpp(1) = 0;

% B0:
B0 = .25*g./t;
B0p = .25*(gp./t-g./t.^2);
B0(1) = .25*(-A/6-.5*B);
B0p(1) = 0;

% A1:
A10 = alph*(A+3*B)/24;
A1 = .125*gp - (1+2*alph)/2*B0 - g.^2/32 - A10;
A1p = .125*gpp - (1+2*alph)/2*B0p - gp.*g/16;
A1p_t = A1p./t;
A1p_t(1) = -A/720 - A^2/576 - A*B/96 - B^2/64 - B/48 + alph*(A/720 + B/48);

% Make f accurately: (Taylor series approx for small t)
fcos = B./(2*cos(t/2)).^2;
f = -A*(1/12 + t.^2/240+t.^4/6048 + t.^6/172800 + t.^8/5322240 + ...
    691*t.^10/118879488000 + t.^12/5748019200 + ...
    3617*t.^14/711374856192000 + 43867*t.^16/300534953951232000);
idx = t > .5;
ti = t(idx);
f(idx) = A.*(1./ti.^2 - 1./(2*sin(ti/2)).^2);
f = f - fcos;

% Integrals for B1: (Note that N isn't large, so we don't need to be fancy).
C = chebcolloc2.cumsummat(N)*(.5*c);
D = chebcolloc2.diffmat(N)*(2/c);
I = (C*A1p_t);
J = (C*(f.*A1));

% B1:
tB1 = -.5*A1p - (.5+alph)*I + .5*J;
tB1(1) = 0;
B1 = tB1./t;
B1(1) = A/720 + A^2/576 + A*B/96 + B^2/64 + B/48 + ...
    alph*(A^2/576 + B^2/64 + A*B/96) - alph^2*(A/720 + B/48);

% A2:
K = C*(f.*tB1);
A2 = .5*(D*tB1) - (.5+alph)*B1 - .5*K;
A2 = A2 - A2(1);

if ( nargout < 3 )
    % Make function for output
    tB1 = @(theta) bary(theta,tB1,t,v);
    A2 = @(theta) bary(theta,A2,t,v);
    return
end

% A2p:
A2p = D*A2;
A2p = A2p - A2p(1);
A2p_t = A2p./t;
% Extrapolate point at t = 0:
w = pi/2-t(2:end);
w(2:2:end) = -w(2:2:end);
w(end) = .5*w(end);
A2p_t(1) = sum(w.*A2p_t(2:end))/sum(w);

% B2:
tB2 = -.5*A2p - (.5+alph)*(C*A2p_t) + .5*C*(f.*A2);
B2 = tB2./t;
% Extrapolate point at t = 0:
B2(1) = sum(w.*B2(2:end))/sum(w);

% A3:
K = C*(f.*tB2);
A3 = .5*(D*tB2) - (.5+alph)*B2 - .5*K;
A3 = A3 - A3(1);

tB1 = @(theta) bary(theta, tB1, t, v);
A2 = @(theta) bary(theta, A2, t, v);
tB2 = @(theta) bary(theta, tB2, t, v);
A3 = @(theta) bary(theta, A3, t, v);
end


function [x, w, v, t] = ultrapts(n, lambda, int, meth)
%ULTRAPTS    Ultraspherical points and Gauss-Gegenbauer quadrature weights.
%   [X]= ULTRAPTS(N, LAMBDA) returns N ultraspherical points X in (-1,1)
%   where the ultraspherical weight function is w(x)=(1-x^2)^(LAMBDA-1/2).
%
%   [X, W] = ULTRAPTS(N, LAMBDA) returns also a row vector W of weights for
%   Gauss-Gegenbauer quadrature.
%   
%   [X, W, V] = ULTRAPTS(N, LAMBDA) returns additionally a column vector V of
%   weights in the barycentric formula corresponding to the points X. The 
%   weights are scaled so that max(abs(V)) = 1.
%
%   ULTRAPTS(N, LAMBDA, INTERVAL) scales the nodes and weights for the 
%   finite interval INTERVAL.
%
%   ULTRAPTS(N, LAMBDA, INTERVAL, METHOD) allows the user to select which 
%   method to use.
%    METHOD = 'REC' uses the recurrence relation for the ultraspherical 
%     polynomials and their derivatives to perform Newton-Raphson iteration. 
%     If 0 < LAMBDA < 1, then the convergence is guaranteed.
%     Default for: LAMBDA <= 3 and N < 100; 3 < LAMBDA <= 8 and N < 500; 
%     8 < LAMBDA <= 13 and N < 1000; 13 < LAMBDA <= 20 and N < 2000;
%     LAMBDA > 20 and N < 3000.
%    METHOD = 'ASY' uses an algorithm adapted from the Hale-Townsend fast 
%     algorithm based upon asymptotic formulae, which is fast and accurate.
%     If 0 < LAMBDA < 1, then the convergence is guaranteed.
%     Default for: LAMBDA <= 3 and N >= 100; 3 < LAMBDA <= 8 and N >= 500; 
%     8 < LAMBDA <= 13 and N >= 1000; 13 < LAMBDA <= 20 and N >= 2000;
%     LAMBDA > 20 and N >= 3000.
%    METHOD = 'GW' will use the traditional Golub-Welsch eigensystem method,
%       which is maintained mostly for historical reasons.
%   
%   [X, W, V, T] = ULTRAPTS(N,LAMBDA) returns also the arccos of the nodes,
%   T = acos(X). In some situations (in particular with 'ASY') these can be
%   computed to a much better relative precision than X.
%
%   The cases LAMBDA=0 and LAMBDA=1 correspond to Gauss-Chebyshev quadratures 
%   nodes and weights, and are treated specially (as a closed form of the nodes 
%   and weights is available). The case LAMBDA=1/2 correspond to Gauss-Legendre
%   quadrature, and it calls LEGPTS which is a more efficient code.
%
% See also CHEBPTS, LEGPTS, JACPTS, LOBPTS, RADAUPTS, HERMPTS, LAGPTS, and
% TRIGPTS.

% Copyright 2017 by The Chebfun Developers. See http://www.chebfun.org/ for
% Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
% 'REC' by Lourenco Peixoto, Jul 2015 - see [3].
% 'ASY' by Lourenco Peixoto, Jan 2016 - see [2] and [4].
%
%  References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature rules",
%       Math. Comp. 23:221-230, 1969.
%   [2] N. Hale and A. Townsend, "Fast and accurate computation of Gauss-Legendre 
%       and Gauss-Jacobi quadrature nodes and weights", SIAM J. Sci. Comp., 2013.
%   [3] L. L. Peixoto, "Desigualdades que garantem a convergencia do metodo
%       de Newton-Raphson para os zeros do polinomio ultraesferico no caso
%       principal", Master's thesis, UFMG, Belo Horizonte, 2015.
%   [4] L. L. Peixoto, "Fast, accurate and convergent computation of 
%       Gauss-Gegenbauer quadrature nodes and weights.", In preparation, 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults:
interval = [-1, 1];
method = 'default';
method_set = 0;

if (lambda <= -0.5)
    error('CHEBFUN:ultrapts:sizeLAMBDA', 'LAMBDA must be greater than -1/2.')
elseif (lambda >= 30)
    warning('CHEBFUN:ultrapts:largeLAMBDA',...
        'LAMBDA >= 30. Results may not be accurate.')
end

% Check inputs:
if ( nargin > 2 )
    if ( nargin == 4 )
        % Calling sequence = ULTRAPTS(N, LAMBDA, INTERVAL, METHOD)
        interval = int;
        method = meth;
        method_set = 1;
    elseif ( nargin == 3 )
        if ( ischar(int) )
            % Calling sequence = ULTRAPTS(N, LAMBDA, METHOD)
            method = int;
            method_set = 1;
        else
            % Calling sequence = ULTRAPTS(N, LAMBDA, INTERVAL)
            interval = int;
        end
    end
    validStrings = {'default', 'GW', 'ASY', 'REC'};
    if ( ~any(strcmpi(method, validStrings)) )
        error('CHEBFUN:ultrapts:inputs', ['Unrecognised input string: ', method]);
    end
    if ( any(isinf(interval)) || numel(interval) ~= 2 || interval(1) >= interval(2) )
        error('CHEBFUN:ultrapts:inputs', 'Interval invalid.');
    end
end


% Deal with trivial cases:
if (n<0)
    error('CHEBFUN:ultrapts:n', 'First input should be positive number.');
elseif (n==0) % Return empty vectors if n==0
    x = [];
    w = [];
    v = [];
    t = [];
    return
elseif (n==1)
    x = 0;
    w = sqrt(pi)*exp(gammaln(lambda+.5)-gammaln(lambda+1));
    v = 1;
    t = 1;
    [x, w] = rescale_ultrapts(x,w,interval,lambda);
    return
elseif (n==2)
    x = [-1; 1]/sqrt(2*(1+lambda));
    w = sqrt(pi)*exp(gammaln(lambda+.5)-gammaln(lambda+1));
    w = [w, w]/2;
    v = [1 ; -1];
    t = acos(x);
    [x, w] = rescale_ultrapts(x,w,interval,lambda);
    return
end

% Special cases:
if ( lambda == 0 ) % Gauss-Chebyshev: lambda = 0
    [x, ~, v] = chebpts(n, interval, 1);
    w = repmat(pi/n,1,n);
    [~, w] = rescale_ultrapts(x, w, interval, lambda);
    t = acos(x);
    return
elseif ( lambda == 1 ) % Gauss-Chebyshev2: lambda = 1
    x = chebpts(n+2, 2);     
    x = x(2:n+1);
    w = pi/(n+1)*(1-x.^2)';
    t = acos(x);
    v = (1-x.^2);  
    v(2:2:end) = -v(2:2:end); 
    [x, w] = rescale_ultrapts(x, w, interval, lambda);
    return
elseif ( lambda == .5 ) % Gauss-Legendre: lambda = 1/2
    [x, w, v, t] = legpts(n, method);
    [x, w] = rescale_ultrapts(x, w, interval, lambda);
    return
end

% Choose the method:
t = [];
if ( lambda <= 3 )
    if ( (n < 100 && ~method_set) || strcmpi(method, 'rec_ultrapts') )
        [x, w, v] = rec_ultrapts(n, lambda); % REC (Recurrence relation) [3]
    elseif ( strcmpi(method, 'GW') )
        [x, w, v] = gw_ultrapts(n,lambda);   % GW [1]
    else
        [x, w, v, t] = asy_ultrapts(n,lambda); % HT [2] and Peixoto [4]
    end
elseif ( lambda <= 8 ) 
    if ( (n < 500 && ~method_set) || strcmpi(method, 'rec_ultrapts') )
        [x, w, v] = rec_ultrapts(n, lambda); % REC (Recurrence relation) [3]
    elseif ( strcmpi(method, 'GW') )
        [x, w, v] = gw_ultrapts(n,lambda);   % GW [1]
    else
        [x, w, v, t] = asy_ultrapts(n,lambda); % HT [2] and Peixoto [4]
    end
elseif ( lambda <= 13 )
    if ( (n < 1000 && ~method_set) || strcmpi(method, 'rec_ultrapts') )
        [x, w, v] = rec_ultrapts(n, lambda); % REC (Recurrence relation) [3]
    elseif ( strcmpi(method, 'GW') )
        [x, w, v] = gw_ultrapts(n,lambda);   % GW [1]
    else
        [x, w, v, t] = asy_ultrapts(n,lambda); % HT [2] and Peixoto [4]
    end
elseif ( lambda <= 20 )
    if ( (n < 2000 && ~method_set) || strcmpi(method, 'rec_ultrapts') )
        [x, w, v] = rec_ultrapts(n, lambda); % REC (Recurrence relation) [3]
    elseif ( strcmpi(method, 'GW') )
        [x, w, v] = gw_ultrapts(n,lambda);   % GW see [1]
    else
        [x, w, v, t] = asy_ultrapts(n,lambda); % HT [2] and Peixoto [4]
    end
else % lambda > 20
    if ( (n < 3000 && ~method_set) || strcmpi(method, 'rec_ultrapts') )
        [x, w, v] = rec_ultrapts(n, lambda); % REC (Recurrence relation) [3]
    elseif ( strcmpi(method, 'GW') )
        [x, w, v] = gw_ultrapts(n,lambda);   % GW see [1]
    else
        [x, w, v, t] = asy_ultrapts(n,lambda); % HT [2] and Peixoto [4]
    end
end
    

% Compute a T is one is asked for:
if ( nargout == 4 && isempty(t) )
    t = acos(x);
end                

% Scale the nodes and quadrature weights:
[x,w] = rescale_ultrapts(x,w,interval,lambda);

% Normalise the barycentric weights:
v = abs(v);
v(2:2:end) = -v(2:2:end); 
v = v./max(abs(v));

end


function [x,w] = rescale_ultrapts(x,w,interval,lambda)
% Rescale to arbitrary finite interval:
if ( ~all(interval == [-1 1]) )
    dab = diff(interval);
    x = (x+1)/2*dab + interval(1);
    w = (.5*dab)^(2*lambda)*w;
end

end
          
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------Routine for GW algorithm--------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v] = gw_ultrapts(n, lambda)
i = (1:n-1)';
bb = .5*sqrt(i.*(i+2*lambda-1)./(i+lambda)./(i+lambda-1));
TT = diag(bb,1) + diag(bb,-1); % Jacobi matrix.
[V, x] = eig( TT ); % Eigenvalue decomposition.
x = diag(x); % Jacobi points.

% Quadrature weights:
w = V(1,:).^2*(2^(2*lambda)*exp(2*gammaln(lambda+.5)-gammaln(2*lambda+1))); 
v = sqrt(1-x.^2).*abs(V(1,:))'; % Barycentric weights.
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------Routine for REC algorithm-------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v] = rec_ultrapts(n,lambda)

% Constants:
lam2 = lambda + lambda;
cte = 2*(n+lambda); % Constant for the weights.
m = floor(.5*(n+1)); % Computes only nonnegative x.
p0 = .5^lambda*exp(.5*gammaln(lam2+1)-gammaln(lambda+.5)); % orthonormal polynomial
% of degree zero.


% Only for nonnegative x.
k = (1:m).';
theta = (k-(1-lambda)*0.5)/(n+lambda)*pi;
x0 = cos(theta);
cos2 = x0.*x0;
% Sharp initial guess (Forster and Petras, 1993):
x = cos( theta + lambda*(1-lambda)/(2*(n+lambda)^2)*(1-(6+lambda*...
    (1-lambda)*(9-2*cos2))./(12*(n+lambda)^2*...
    (1-cos2))).*cot(theta) );

if (lambda>0 && lambda<1)
    % Choose initial guess for guaranteed convergence [3]:
    if (x > x0)
        x = x0;
    end
elseif ( n > 21 )
    % Gatteschi's approximation (1979) for initial guesses on the
    % boundary region:
    Lam = lambda*(1-lambda);
    N = sqrt((n+lambda)^2+lambda*(1-lambda)/3);
    bz = besselroots(lambda-.5,10); % Approximates zeros of Bessel.
    x(1:10) = cos(bz/N-Lam/90.*(bz.^3+2*(lambda^2-lambda-.75).*bz)/N^5);
end

dx = inf;
counter = 0;

r1=zeros(n,1);
r2=r1;

% The terms of recurrence relation for the orthonormal polynomials:
for k = 1:n
    r1(k) = 2*sqrt((k-1+lambda)*(k+lambda)/(k*(k-1+lam2)));
    r2(k) = sqrt((k-1)*(k+lambda)*(k+lam2-2)/(k*(k+lam2-1)*(k+lambda-2)));
end

% Loop until convergence:
while ( norm(dx, inf) > eps && counter < 20 )
    % Initialise:
    P = p0;
    P1 = 0;
    counter = counter + 1;
    for k = 1:n,
        P2 = P1;
        P1 = P;
        P = r1(k)*x.*P1 - r2(k)*P2; % P(x) is the orthonormal polynomial.
    end
    PP = ( sqrt((n+lam2-1)*n*(n+lambda)/(n-1+lambda))*P1 - n*x.*P )...
        ./(1-x.*x);
    % Newton step:
    dx = P./PP;
    % Newton update:
    x = x - dx;
end

P = p0;
P1 = 0;

 % Once more for P1 and PP:
 for k = 1:n,
     P2 = P1;
     P1 = P;        
     P = r1(k)*x.*P1 - r2(k)*P2;
 end
 c2 = 1-x.*x;
 PP = ( sqrt((n+lam2-1)*n*(n+lambda)/(n-1+lambda))*P1 - n*x.*P )...
         ./(c2);

% Quadrature weights for nonnegative values:
 w = cte./(c2.*PP.*PP); % Usual relation. 

% Reflect for negative values:
s = mod(n,2);

x = [-x(1:m-s); x(m:-1:1)];

% Reflect for quadrature weights:
w = [w(1:m-s) ; w(m:-1:1)]';

% Reflect for derivatives:
ders = [PP(1:m-s); PP(m:-1:1)];

% Barycentric weights:
v = 1./ders;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------Routine for ASY algorithm------------------------ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w,v,t] = asy_ultrapts(n,lambda)
% ASY computes nodes and weights using asymptotic formulae.

if (n<=21) % Use only interior formula:
    nbdy = floor(.5*(n+1));
    [x, w, v, t] = asy1_ultrapts(n, lambda, nbdy);
    return
end

% Determine switch between interior and boundary regions:
nbdy = 10; % Typically, the 10 nodes nearest the boundary.

% Interior algorithm:
[x, w, v, t] = asy1_ultrapts(n, lambda, nbdy);

% Boundary algorithm:
[x2, w2, v2, t2] = asy2_ultrapts(n, lambda, nbdy);

% Combine:
bdyidx1 = n-(nbdy-1):n;
bdyidx2 = nbdy:-1:1;
x(bdyidx1) = x2;
w(bdyidx1) = w2;
v(bdyidx1) = v2;
t(bdyidx1) = t2;

% Reflect using symmetry:
x(bdyidx2) = -x2;
w(bdyidx2) = w2;
v(bdyidx2) = v2;
t(bdyidx2) = pi-t2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              ASY1 (interior)                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w,v,t] = asy1_ultrapts(n,lambda,nbdy)
% Algorithm for computing nodes and weights in the interior region.
% N > 21.

k = (floor((n+1)/2):-1:1);
 
t0 = (k-(1-lambda)*0.5)/(n+lambda)*pi;
cos2 = cos(t0).^2;
% Sharp initial guess (Forster and Petras, 1993):
t = t0 + lambda*(1-lambda)/(2*(n+lambda)^2)*(1-(6+lambda*...
    (1-lambda)*(9-2*cos2))./(12*(n+lambda)^2*...
    (1-cos2))).*cot(t0);
if (lambda>0 && lambda<1)
    % Choose initial guess for guaranteed convergence [4]:
    if (t < t0)
        t = t0;
    end
end


% Locate the boundary node:
mint = t(end-nbdy+1);
idx = max(find(t<mint,1)-1,1);

% Initialise:
dt = inf;
% Newton iteration:
while ( norm(dt, inf) > sqrt(eps)/1000)       % <-- Enough, as again below
    [vals, ders] = feval_asy1_ultrapts(n, t, 1);        % Evaluate via asy_ultrapts formulae
    dt = vals./ders;                           % Newton step
    t = t - dt;
    dt = dt(1:idx-1);                          % Ignore boundary terms
end
[vals, ders] = feval_asy1_ultrapts(n, t, 1);  % Once more for good ders. 
t = transpose(t - vals./ders);


% Constant for the weights: 
% ( cte = 4^(-lambda)*pi*gamma(n+2*lambda)*gamma(n+1)/gamma(n+lambda)^2 )
% See ratio of gamma functions in [2].
dt = .5*(lambda-1)^2/n;
ds = .5*lambda^2/(n+lambda-1);
s = ds - dt;
j = 1;
while ( abs((ds-dt)/s) > eps/100 ) % Taylor series in expansion
    j = j+1;
    k = (j-1)/(j+1);
    dt = -(lambda-1)*k/n*dt;
    ds = -lambda/(n+lambda-1)*k*ds;
    s = s+ds-dt;
end

p2 = exp(s)*sqrt(n*(n+2*lambda-1))*(1+(lambda-1)/n)^(lambda-1);
% Stirling's series:
g = [1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880, ...
     5246819/75246796800, -534703531/902961561600, ...
     -4483131259/86684309913600, 432261921612371/514904800886784000];
f = @(z) sum(g.*[1, cumprod(ones(1, 9)./z)]);
cte = 4^(-lambda)*pi*p2*f(n+2*lambda-1)*f(n)/(f(n+lambda-1)^2);

% Compute x, w, and v:
x = cos(t);
w = (cte./ders.^2);
v = (transpose(sin(t))./ders).';

if (mod(n,2))
    x(1) = 0; % computed analytically
    t(1) = pi/2; % computed analytically
    x = [-x(end:-1:2); x];
    w = [w(end:-1:2), w];
    v = -[v(end:-1:2); v];
    t = [pi-t(end:-1:2); t];
else
    x = [-x(end:-1:1); x];
    w = [w(end:-1:1), w];
    v = [-v(end:-1:1); v];
    t = [pi-t(end:-1:1); t];
end


function [vals, ders] = feval_asy1_ultrapts(n, t, flag) 
% Evaluate asymptotic formula (interior) - Szego p. 197 (8.21.14). 

% The polynomial is normalised in order to avoid the computation of 
% gamma(lambda) into the weights:
% If P is the ultraspherical polynomial, and P_norm is the normalised, then
% P_norm = ( gamma(n+1)*gamma(lambda)/(2*gamma(n+lambda)) ) * P.

% Number of expansion terms:
r = mod(lambda,2);
if r==0 || r==1
    % If lambda is an integer, then M=lambda-1 gives accurate formula! 
    M = lambda-1;
else
    % M = 30 on otherwise. (Obtained experimentally.)
    M = 30;
end

% Coefficients in expansion:
c = cumprod( (1-lambda:M-lambda)./(n+lambda-1:-1:n+lambda-M) );
d = cumprod((1+lambda:M-1+lambda)./(2:M));
d = lambda*[1, d];
c = [1, c.*d];

% Total number of expansion terms is M+1: 
M=M+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Some often used vectors/matrices:
onesT = ones(1, length(t));
onesM = ones(M, 1);
Mlam = transpose((0:M-1)+lambda);
Mnlam = transpose(n+lambda-(0:M-1));
onesMcotT = onesM*cot(t);
MnlamonesT = Mnlam*onesT;
MlamonesT = Mlam*onesT;

twoSinT = onesM*(2*sin(t));
denom = cumprod(twoSinT)./(twoSinT).^(1-lambda);

if ( ~flag )
    alpha = MnlamonesT.*(onesM*t) - .5*pi*MlamonesT;
    cosAlpha = cos(alpha);
    sinAlpha = sin(alpha);
else
    % Adapted from JACPTS, and LEGPTS:
    %%%%%%%%%%%%%%%% Taylor expansion of cos(alpha0) %%%%%%%%%%%%%%
    k = numel(t):-1:1;
    % HI-LO expansion, to accurately compute (n+lambda)*t - (k-.5*lambda)*pi
    ta = double(single(t));
    tb = t - ta;
    hi = n*ta;
    lo = n*tb+lambda*t;
    pia = double(single(pi)); 
    pib = pi - pia;
    dh = (hi-(k-.25)*pia) + lo -.5*(lambda-.5)*pia -...
        (k-.25+.5*(lambda-.5))*pib;

    % Compute cosAlpha(1,:) using Taylor series:
    tmp = 0; sgn = 1; fact = 1; DH = dh; dh2 = dh.*dh;
    for jj = 0:100
        dc = sgn*DH/fact;
        tmp = tmp + dc;
        sgn = -sgn;
        fact = fact*(2*jj+3)*(2*jj+2);
        DH = DH.*dh2;
        if ( norm(dc,inf) ) < eps/2000, break, end
    end
    tmp(2:2:end) = -tmp(2:2:end);
    [~, loc] = max(abs(tmp));
    tmp = sign(cos((n+lambda)*t(loc)-.5*lambda*pi)*tmp(loc))*tmp;
    cosAlpha(1,:) = tmp;

    % Compute sinAlpha(1,:) using Taylor series:
    tmp = 0; sgn = 1; fact = 1; DH = 1; dh2 = dh.*dh;
    for jj = 0:100
        dc = sgn*DH/fact;
        tmp = tmp + dc;
        sgn = -sgn;
        fact = fact*(2*jj+2)*(2*jj+1);
        DH = DH.*dh2;
        if (norm(dc, inf)) < eps/2000, break, end
    end
    tmp(2:2:end) = -tmp(2:2:end);
    [~, loc] = max(abs(tmp));
    tmp = sign(sin((n+lambda)*t(loc)-.5*lambda*pi)*tmp(loc))*tmp;
    sinAlpha(1,:) = tmp;

    % Compute cosAlpha(k,:) and sinAlpha(k,:) for k = 2,...,M:
    sint = sin(t);
    cost = cos(t);
    for kk = 2:M
        cosAlpha(kk,:) = sinAlpha(kk-1,:).*cost - cosAlpha(kk-1,:).*sint;
        sinAlpha(kk,:) = -(cosAlpha(kk-1,:).*cost + sinAlpha(kk-1,:).*sint);
    end
end

% Sum up all the terms:
vals = c*(cosAlpha./denom); % P(theta)
numer = MnlamonesT.*sinAlpha + MlamonesT.*cosAlpha.*onesMcotT;
ders = -c*(numer./denom); % (dP/dtheta)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              ASY2 (boundary)                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v, t] = asy2_ultrapts(n, lambda, nbdy) 
% Algorithm for computing nodes and weights near the boundary.
% N > 2.

k=(1:nbdy).';

% Constant:
Lam = lambda*(1-lambda);

% Choose initial guess:
if ( lambda > 1 || lambda < 0 )
    % Gatteschi's approximation (1979):
    N = sqrt((n+lambda)^2+lambda*(1-lambda)/3);
    bz = besselroots(lambda-.5,nbdy); % Approximates zeros of Bessel
    t = bz/N-Lam/90.*(bz.^3+2*(lambda^2-lambda-.75).*bz)/N^5;
else
    % This initial guess guarantees convergence for 0<lambda<1 [4]:
    t0 = (k-(1-lambda)*0.5)/(n+lambda)*pi;
    cos2 = cos(t0).^2;
    % Sharp initial guess (Forster and Petras, 1993):
    t = t0 + lambda*(1-lambda)/(2*(n+lambda)^2)*(1-(6+lambda*...
        (1-lambda)*(9-2*cos2))./(12*(n+lambda)^2*...
        (1-cos2))).*cot(t0);
    % Choose initial guess for guaranteed convergence [4]:
    if (t < t0)
        t = t0;
    end
end

% Useful constants:
a=lambda-.5;
b=a;

dt = inf; j = 0;
% Newton iteration: 
while ( norm(dt,inf) > sqrt(eps)/1000 && j < 20) 
    [tB1, A2, tB2, A3, Ak, Cint, Dint, f2, v] = asy2_firstterms_ultrapts(a, b, t, n);
    [vals, ders] = feval_asy2_ultrapts(n, t, 0); % Evaluate via asymptotic formula.
    dt = vals./ders;                    % Newton update.
    t = t + dt   ;                      % Next iterate.
    j = j + 1;
end
[tB1, A2, tB2, A3, Ak, Cint, Dint, f2, v] = asy2_firstterms_ultrapts(a, b, t, n);
[vals, ders] = feval_asy2_ultrapts(n, t, 1);     % Evaluate via asymptotic formula.
dt = vals./ders;                        % Newton update
t = t + dt;    
    
% flip:
t = t(nbdy:-1:1); 
ders = ders(nbdy:-1:1);

% Constant for the weights:
% (cte = (2^(2*lambda+1)*(n+lambda)^(2*lambda-1)*gamma(n+1)/gamma(n+2*lambda)
% See ratio of gamma functions in [2].
ds = -.5*(2*lambda-1)^2/n;
s = ds;
j = 1;
while ( abs(ds/s) > eps/1000 ) % Taylor series in expansion 
    j = j+1;
    ds = -(2*lambda-1)*(j-1)/(j+1)/n*ds;
    s = s + ds;
end
p2 = exp(s)/sqrt(n+2*lambda-1);
% Stirling's series:
g = [1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880, ...
    5246819/75246796800, -534703531/902961561600, ...
    -4483131259/86684309913600, 432261921612371/514904800886784000];
f = @(z) sum(g.*[1, cumprod(ones(1, 9)./z)]);
cte = 2^(2*lambda+1)*(1+lambda/n)^(2*lambda-1)*sqrt(n)*p2*(f(n)/f(n+2*lambda-1));

% Revert to x-space:
x = cos(t);      
w = (cte./ders.^2).';
if ( find(w==0) )
    warning('CHEBFUN:ultrapts:largeNLAMDBA',...
       'Some WEIGHTS near the boundary region become zero due to overflow.');
end % The overflow occurs on the computation of ders.^2  


% Constant for the barycentric weights:
% C1 = gamma(n+2*lambda)/gamma(n+lambda)
% See ratio of gamma functions in [2].
m = n+lambda;
ds = .5*lambda^2/(m-1);
s = ds;
j = 1;
while ( abs(ds/s) > eps/100 ) % Taylor series in expansion 
    j = j+1;
    ds = -lambda*(j-1)/(j+1)/(m-1)*ds;
    s = s + ds;
end
p2 = exp(-s)*(n+2*lambda-1)^(-.5)*(1-1/(n+lambda))^(.5-lambda);
% Stirling's series:
g = [1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880, ...
    5246819/75246796800, -534703531/902961561600, ...
    -4483131259/86684309913600, 432261921612371/514904800886784000];
f = @(z) sum(g.*[1, cumprod(ones(1, 9)./z)]);
C1 = 2^(2*lambda+.5)/sqrt(pi)*p2*(f(m-1)/f(m+lambda-1));

% Revert to x-space:
v = sin(t)./(ders/C1); % barycentric weights with dP/dtheta computed such 
% as in interior region.


% The function below is adapted from JACPTS:
function [vals, ders] = feval_asy2_ultrapts(n, t, flag)
% Evaluate asymptotic formula (boundary) - Baratella and Gatteschi (1998).
% It also computes additional terms for asymptotic series as needed. 

% The polynomial is normalised in order to avoid the computation of 
% gamma(lambda) on the weights:
% If P is the ultraspherical polynomial, and P_norm is the normalised, then
% P_norm = ( sqrt(2)*(n+lambda)^(lambda-.5)*gamma(2*lambda)*gamma(n+1) / 
% gamma(lambda+.5)*gamma(n+2*lambda) ) * P.
    
% Useful constants:
rho = n + lambda;
rho2 = rho - 1;
A = Lam;
        
% Evaluate the Bessel functions:
Ja = besselj(a, rho*t, 0);
Jb = besselj(a + 1, rho*t, 0);
Jbb = besselj(a + 1, rho2*t, 0);
if ( ~flag )
    Jab = besselj(a, rho2*t, 0);
else
    % In the final step, perform accurate evaluation
    Jab = besselTaylor_ultrapts(-t, rho*t, a);
end
% Evaluate functions for recursive definition of coefficients:
gt = 2*A*(cot(t)-1./t);
gtdx = .5*A*(4./t.^2-csc(t/2).^2-sec(t/2).^2);
tB0 = .25*gt;
A10 = a*A/6;
A1 = gtdx/8 - (1+2*a)/8*gt./t - gt.^2/32 - A10;

% Higher terms:
tB1t = tB1(t); 
A2t = A2(t); 

% VALS:
vals = Ja + Jb.*tB0/rho + Ja.*A1/rho^2 + Jb.*tB1t/rho^3 + Ja.*A2t/rho^4;
% DERS:
ders = Jab + Jbb.*tB0/rho2 + Jab.*A1/rho2^2 + Jbb.*tB1t/rho2^3 +...
    Jab.*A2t/rho2^4;

% Higher terms:
tB2t = tB2(t); A3t = A3(t);
dv = Jb.*tB2t/rho^5 + Ja.*A3t/rho^6;
vals = vals + dv;
dd = Jbb.*tB2t/rho2^5 + Jab.*A3t/rho2^6;
ders = ders + dd;
   
% Scaling:
valstmp = vals;
denom = (sin(t)/2).^lambda;
vals = (sqrt(t)./denom).*valstmp; % P(theta)

% Relation for derivative:
C2 = n/(n+a)*(rho/rho2)^a;
ders = (-n*(2*(n+lambda)-1)*cos(t).*valstmp + 2*(n+a)^2*C2*ders)/...
    (2*(n+lambda)-1);
ders = ders.*(sqrt(t)./(denom.*sin(t))); % dP(theta)

% Increase terms as needed:
if ( lambda > 1 || lambda < 0 )
    
    % Initialise:
    del = inf;
    deld = del;
    
    k = 3; % Ak = A3 - A3(1)
    alph = lambda-.5; % constant
    
    while ( del > eps && deld > eps && k<=60 )
        % Akp:
        Akp = Dint*Ak;
        Akp = Akp - Akp(1);
        Akp_t = Akp./t;
        % Extrapolate point at t = 0:
        w = pi/2-t(2:end);
        w(2:2:end) = -w(2:2:end);
        w(end) = .5*w(end); 
        Akp_t(1) = sum(w.*Akp_t(2:end))/sum(w);

        % Bk:
        tBk = -.5*Akp - (.5+alph)*(Cint*Akp_t) + .5*Cint*(f2.*Ak);
        Bk = tBk./t;
        % Extrapolate point at t = 0
        Bk(1) = sum(w.*Bk(2:end))/sum(w);

        % Ak1:
        K = Cint*(f2.*tBk);
        Ak1 = .5*(Dint*tBk) - (.5+alph)*Bk - .5*K;
        Ak1 = Ak1 - Ak1(1);
        Atemp = Ak1;

        tBk = @(theta) bary(theta,tBk,t,v);
        Ak1 = @(theta) bary(theta,Ak1,t,v);

        tBkt = tBk(t); Ak1t = Ak1(t);

        dv = Jb.*tBkt/rho^(2*k+1) + Ja.*Ak1t/rho^(2*k+2);
        dd = Jbb.*tBkt/rho2^(2*k+1) + Jab.*Ak1t/rho2^(2*k+2);

        del = sqrt(t).*dv./denom;
        vals = vals + del; % P(theta)

        deld = ((-n*(2*(n+lambda)-1)*cos(t).*dv + 2*(n+a)^2*C2*dd)...
            /(2*(n+lambda)-1)).*(sqrt(t)./(denom.*sin(t)));
        ders = ders + deld; % dP(theta)

        k = k + 1;
        Ak = Atemp;
        del = norm(del./vals,inf);
        deld = norm(deld./ders,inf);
    end
    
end

end

end

% This code is from JACPTS:
function Ja = besselTaylor_ultrapts(t, z, a)
% BESSELTAYLOR    Accurate evaluation of Bessel function J_A for asy2_ultrapts. (See [2].)
% BESSELTAYLOR(T, Z, A) evaluates J_A(Z+T) by a Taylor series expansion about Z. 

npts = numel(t);
kmax = 30;
H = bsxfun(@power, t, 0:kmax).';
% Compute coeffs in Taylor expansions about z (See NIST 10.6.7)
[nu, JK] = meshgrid(-kmax:kmax, z);
[Bjk] = besselj(a + nu, JK, 0);
nck = abs(pascal(floor(1.25*kmax), 1)); nck(1,:) = []; % nchoosek
AA = [Bjk(:,kmax+1), zeros(npts, kmax)];
fact = 1;
for k = 1:kmax
    sgn = 1;
    for l = 0:k
        AA(:,k+1) = AA(:,k+1) + sgn*nck(k,l+1)*Bjk(:,kmax+2*l-k+1);
        sgn = -sgn;
    end
    fact = k*fact;
    AA(:,k+1) = AA(:,k+1)/2^k/fact;
end
% Evaluate Taylor series:
Ja = zeros(npts, 1);
for k = 1:npts
    Ja(k,1) = AA(k,:)*H(:,k);
end
end


%This code is adapted from JACPTS:
function [tB1,A2,tB2,A3,Ak,C,D,f,v] = asy2_firstterms_ultrapts(alph, bet, theta, n)
% ASY2_FIRSTTERMS   First terms for boundary asymptotic series.
% Compute the first terms in asy2_ultrapts boundary formula.

% These constants are more useful than alph and bet:
A = (0.25 - alph^2);
B = (0.25 - bet^2);

% For now, just work on half of the domain:
 c = max(max(theta), .5);
if ( n < 30 )
    N = ceil(40 - n);
elseif ( n >= 30 && c > pi/2-.5)
    N = 15;
else
    N = 10;
end
Nm1 = N - 1;


% Scaled 2nd-kind Chebyshev points and barycentric weights:
t = .5*c*( sin(pi*(-Nm1:2:Nm1)/(2*Nm1)).' + 1 );
v = [.5 ; ones(Nm1,1)];
v(2:2:end) = -1;
v(end) = .5*v(end);

% The g's:
g = A*(cot(t/2)  -2./t) - B*tan(t/2);
gp = A*(2./t.^2 - .5*csc(t/2).^2) - .5*(.25-bet^2)*sec(t/2).^2;
gpp = A*(-4./t.^3 + .25*sin(t).*csc(t/2).^4) - 4*B*sin(t/2).^4.*csc(t).^3;
g(1) = 0; gp(1) = -A/6-.5*B; gpp(1) = 0;

% B0:
B0 = .25*g./t;
B0p = .25*(gp./t-g./t.^2);
B0(1) = .25*(-A/6-.5*B);
B0p(1) = 0;

% A1:
A10 = alph*(A+3*B)/24;
A1 = .125*gp - (1+2*alph)/2*B0 - g.^2/32 - A10;
A1p = .125*gpp - (1+2*alph)/2*B0p - gp.*g/16;
A1p_t = A1p./t;
A1p_t(1) = -A/720 - A^2/576 - A*B/96 - B^2/64 - B/48 + alph*(A/720 + B/48);

% Make f accurately: (Taylor series approx for small t)
fcos = B./(2*cos(t/2)).^2;
f = -A*(1/12 + t.^2/240+t.^4/6048 + t.^6/172800 + t.^8/5322240 + ...
    691*t.^10/118879488000 + t.^12/5748019200 + ...
    3617*t.^14/711374856192000 + 43867*t.^16/300534953951232000);
idx = t > .5;
ti = t(idx);
f(idx) = A.*(1./ti.^2 - 1./(2*sin(ti/2)).^2);
f = f - fcos;

% Integrals for B1: (Note that N isn't large, so we don't need to be fancy).
C = chebcolloc2.cumsummat(N)*(.5*c);
D = chebcolloc2.diffmat(N)*(2/c);
I = (C*A1p_t);
J = (C*(f.*A1));

% B1:
tB1 = -.5*A1p - (.5+alph)*I + .5*J;
tB1(1) = 0;
B1 = tB1./t;
B1(1) = A/720 + A^2/576 + A*B/96 + B^2/64 + B/48 + ...
    alph*(A^2/576 + B^2/64 + A*B/96) - alph^2*(A/720 + B/48);

% A2:
K = C*(f.*tB1);
A2 = .5*(D*tB1) - (.5+alph)*B1 - .5*K;
A2 = A2 - A2(1);

% A2p:
A2p = D*A2;
A2p = A2p - A2p(1);
A2p_t = A2p./t;
% Extrapolate point at t = 0:
w = pi/2-t(2:end);
w(2:2:end) = -w(2:2:end);
w(end) = .5*w(end);
A2p_t(1) = sum(w.*A2p_t(2:end))/sum(w);

% B2:
tB2 = -.5*A2p - (.5+alph)*(C*A2p_t) + .5*C*(f.*A2);
B2 = tB2./t;
% Extrapolate point at t = 0:
B2(1) = sum(w.*B2(2:end))/sum(w);

% A3:
K = C*(f.*tB2);
A3 = .5*(D*tB2) - (.5+alph)*B2 - .5*K;
A3 = A3 - A3(1);

% temporary term:
Ak = A3;

% Make function for output:
tB1 = @(theta) bary(theta, tB1, t, v);
A2 = @(theta) bary(theta, A2, t, v);
tB2 = @(theta) bary(theta, tB2, t, v);
A3 = @(theta) bary(theta, A3, t, v);
end
