function varargout = pwlinstab(f1, f2, phi, x, Q, opts)
% Performs a linear stability analysis for a piecewise polynomial system,
%
%          / f1(x)      if phi(x) <= 0;
%   xdot = |
%          \ f2(x)      else;
%
% about the equilibrium point x=0. The piecewise linearizations are 
% respectively given as
% 
%   xdot = Ai*x. 
%
%% Usage & description
%
%   [V,A1,A2,P] = pwlinstab(f1,f2,phi,x)
%   [V1,V2,...] = pwlinstab(f1,f2,phi,x)
%   [...] = pwlinstab(...,Q)
%   [...] = pwlinstab(...,Q,opts)
%
% Inputs:
%       -f1:  first polynomial vector field
%       -f2:  second polynomial vector field
%       -phi: boundary condition (scalar field)
%       -x:   state-space vector as PVAR
%       -Q:   postive definite parameter in the Lyapunov equations; must be
%             either scalar or square matrix of the size of |x|, or empty. 
%             [default = 1e-6]
%       -opts:  SOS options structre; see SOSOPTIONS.
%
% Outputs:
%       -V:   quadratic Lyapunov function, V = x'*P*x
%       -A1,A2:  linearizations of the piecewise systems around x=0
%       -P:   positive definite solution of the Lyapunov equations; this 
%             matrix is used to form the quadratic Lyapunov function. 
%
% If the linearization is not stable then V and P are returned as empty.
%
%%

assert(nargout ~= 3 && nargout <= 5, 'Undefined number of outputs (%g)', nargout);


if ~exist('Q', 'var') || isempty(Q)
    Q = 1e-6;
end
if ~exist('opts', 'var')
    opts = sosoptions;
end

varargout = cell(1,nargout);
P = [];

% Linearize: xdot = A*x
A1 = plinearize(f1,x);
A2 = plinearize(f2,x);

ev1 = eig(A1);
ev2 = eig(A2);
if max(real(ev1), real(ev2)) >= 0
    % if A1, A2 are unstable
    % nothing to do
elseif nargout <= 1 || nargout == 4
    % if A1, A2 are stable solve LMI
    %   A1'*P + P*A1 + S1 < 0
    %   A2'*P + P*A2 + S2 < 0
    % for common Lyapunov function
    [V,P] = sosdecvar('p', x);
    
    sosc = polyconstr;
    sosc(1) = x'*(A1'*P + P*A1)*x - phi <= -x'*Q*x;
    sosc(2) = x'*(A2'*P + P*A2)*x + phi <= -x'*Q*x;
    sosc(3) = V >= x'*Q*x;
    
    [info,dopt] = sosopt(sosc, x, opts);

    if info.feas
        P = subs(P, dopt);
        varargout{1} = x'*P*x;
    end
else
    % if A1, A2 are stable solve LMI
    %   A1'*P1 + P1*A1 + S1 < 0
    %   A2'*P2 + P2*A2 + S2 < 0
    % for multiple Lyapunov function
    [V1,P1] = sosdecvar('p1', x);
    [V2,P2] = sosdecvar('p2', x);
    
    % continuity decision variables
    r1 = sosdecvar('c1', x);
    r2 = sosdecvar('c2', x);
    r3 = sosdecvar('c3', x);
    r4 = sosdecvar('c4', x);
    
    sosc = polyconstr;
    sosc(1) = x'*(A1'*P1 + P1*A1)*x - phi <= -x'*Q*x;
    sosc(2) = x'*(A2'*P2 + P2*A2)*x + phi <= -x'*Q*x;
    sosc(3) = V1 >= x'*Q*x;
    sosc(4) = V2 >= x'*Q*x;
    
    sosc(5) = -((V1-V2) - r1*phi + r2*phi) >= 0;
    sosc(6) = -((V2-V1) - r3*phi + r4*phi) >= 0;
    
    sosc(7)  = r1 >= 0;
    sosc(8)  = r2 >= 0;
    sosc(9)  = r3 >= 0;
    sosc(10) = r4 >= 0;
    
    [info,dopt] = sosopt(sosc, x, opts);
    
    if info.feas
        P = {
            subs(P1,dopt)
            subs(P2,dopt)
        };
        varargout{1} = x'*P{1}*x;
        varargout{2} = x'*P{2}*x;
    end
end

if nargout >= 4
    varargout{end-2} = A1;
    varargout{end-1} = A2;
    varargout{end}   = P;
end
