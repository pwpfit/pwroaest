function [V, A1, A2, P] = pwlinstab(f1, f2, phi, x, Q, opts)
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

if ~exist('Q', 'var') || isempty(Q)
    Q = 1e-6;
end
if ~exist('opts', 'var')
    opts = sosoptions;
end

% Linearize: xdot = A*x
A1 = plinearize(f1,x);
A2 = plinearize(f2,x);

ev1 = eig(A1);
ev2 = eig(A2);
if min(real(ev1), real(ev2)) < 0
    % if A1, A2 are stable solve LMI
    %   A1'*P + P*A1 = -Q
    %   A2'*P + P*A2 = -Q
    [V,P] = sosdecvar('p', x);
    
    sosc = polyconstr;
    sosc(1) = x'*(A1'*P + P*A1)*x - phi <= -x'*Q*x;
    sosc(2) = x'*(A2'*P + P*A2)*x + phi <= -x'*Q*x;
    sosc(3) = V >= x'*Q*x;
    
    [info,dopt] = sosopt(sosc, x, opts);
else
    info.feas = false;
end

if info.feas
    P = subs(P, dopt);
    V = x'*P*x;
else
    V = [];
    P = [];
end