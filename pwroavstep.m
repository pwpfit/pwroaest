function [V,c] = pwroavstep(f1,f2,phi,p,x,z,beta,gamma,s0,s,si,varargin)
% Solves the V-s of the piecewise ROA iteration.
%
%% Usage & description
%
%   [V,c] = pwroavstep(f1,f2,phi,p,x,z,beta,gamma,s0,s2,si)
%   [...] = pwroavstep(...,L1,L2,opts)
%
% Solves for a common Lyapunov function V which satisfies the SOS
% constraints while holding all other variables fixed.
%
% Inputs:
%       -f1:  first polynomial vector field
%       -f2:  second polynomial vector field
%       -phi: boundary condition (scalar field)
%       -p:   shape function (scalar field)
%       -x:   state-space vector as PVAR
%       -z:   Nz-by-1 column vector of monomials; specifies the Lyapunov
%             function decision variable in the vector form V(x) = c'*z(x).
%       -beta:  level set of shape function
%       -gamma: level set of Lyapunov function
%       -s1:  multiplier for beta-step 
%       -s2:  multipliers [s2'; s2"] for gamma-step (Lyapunov gradient)
%       -si:  multipliers [si'; si"] for gamma-step (boundary condition)
%       -L1:  epsilon 1 (double or scalar field); enforces strict positive 
%             definiteness of the Lyapunov function.
%       -L2:  epsilon 2 (double or scalar field); enforces strict negative
%             definiteness of the Lyapunov gradients.
%       -opts:  SOS options structre; see SOSOPTIONS.
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
% * Created:    2018-05-23
% * Changed:    2018-05-23
%
%% See also
%
% See PWROAEST, ROAVSTEP
%%

for i=1:length(varargin)
    arg = varargin{i};
    
    if ~exist('L1', 'var'),       L1 = arg;
    elseif ~exist('L2', 'var'),   L2 = arg;
    elseif ~exist('opts', 'var'), opts = arg;
    else
        error('Undefined input for pwroavstep: %s.', class(arg));
    end
end
if ~exist('L1','var'), L1 = 0; end
if ~exist('L2','var'), L2 = 0; end
if ~exist('opts','var')
    opts = sosoptions;
end

% Lyapunov decision variable
[V,c] = polydecvar('c',z);

%% V-s feasibility problem
sosconstr = polyconstr;

% V-L1 in SOS
sosconstr(1) = V >= L1;

% {x: p(x) <= b} is contained in {x: V(x) <= g}
sosconstr(2) = -((V-gamma) + s0*(beta-p)) >= 0;

% {x: V(x) <= g} is contained in {x: grad(V)*f < 0}
gradV = jacobian(V,x);
% -( pa + (g-p2)*s - phi*si ) in SOS
sosconstr(3) = -(gradV*f1 + L2 + s(1)*(gamma-V) - si(1)*phi) >= 0;
% -( pb + (g-p2)*s + (phi-l)*si ) in SOS
sosconstr(4) = -(gradV*f2 + L2 + s(2)*(gamma-V) + si(2)*(phi-L2)) >= 0;

% solve problem
[info,dopt] = sosopt(sosconstr,x,opts);

%% Output
if info.feas
    V = subs(V,dopt);
    c = subs(c,dopt);
else
    V = [];
    c = [];
end
