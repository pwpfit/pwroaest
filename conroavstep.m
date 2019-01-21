function [V,c] = conroavstep(f,con,p,x,z,beta,gamma,s1,s2,sg,L1,L2,opts)
% Solves the V-s step of the piecewise ROA iteration.
%
%% Usage & description
%
%   [V,c] = conroavstep(f,c,p,x,z,beta,gamma,s1,s2,sg)
%   [...] = conroavstep(...,L1,L2,opts)
%
% Solves for a Lyapunov function V which satisfies the SOS constraints 
% while holding all other variables fixed.
%
% Inputs:
%       -f:   polynomial vector field
%       -con: constraint function (scalar polynomial)
%       -p:   shape function (scalar polynomial)
%       -x:   state-space vector as PVAR
%       -z:   Nz-by-1 column vector of monomials; specifies the Lyapunov
%             function decision variable in the vector form V(x) = c'*z(x).
%       -beta:     level set of shape function
%       -gamma:    level set of Lyapunov function
%       -s1,s2,sg: set containment multipliers
%       -L1:  epsilon 1 (double or scalar field); enforces strict positive 
%             definiteness of the Lyapunov function.
%       -L2:  epsilon 2 (double or scalar field); enforces strict negative 
%             definiteness of the Lyapunov gradient.
%       -opts:     SOS options structure;  see SOSOPTIONS.
%
% Outputs:
%       -V: Lyapunov function satisfying the SOS constraints; will be empty
%           if no feasible solution is found.
%       -c: Nz-by-1 vector of coefficients of the returned Lyapunov
%           function.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2019-01-21
% * Changed:    2019-01-21
%
%% See also
%
% See PWROAEST, ROAVSTEP
%%

if isempty(opts)
    opts = sosoptions;
end


% Create Lyapunov decision variable
[V,c] = polydecvar('c',z);

% Constraint 1: 
% V-L1  in SOS
sosconstr = polyconstr;
sosconstr(1) = V >= L1;

% Constraint 2: 
% {x: p(x) <= b} is contained in {x: V(x) <= g}
sosconstr(2) = (V-gamma) <= s1*(p-beta);

% Constraint 3:
% {x: V(x) <= g} is contained in {x: grad(V)*f < 0}
Vdot = jacobian(V,x)*f;
sosconstr(3) = Vdot <= -L2 + s2*(V-gamma);

% Constraint 4:
% {x: V(x) <= g} is contained in {x: c(x) <= 0}
sosconstr(4) = -(con + sg*(gamma-V)) >= 0;

% Solve feasibility problem
[info,dopt] = sosopt(sosconstr,x,opts);

% Create output
if info.feas
    V = subs(V,dopt);
    c = subs(c,dopt);
else
    V=[];
    c=[];
end

