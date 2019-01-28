function [K,c] = roaKstep(f,con,~,x,u,z,V,~,gamma,~,s2,sg,~,L2,opts)
% Solves the K-V-s step of the constraint ROA iteration.
%
%% Usage & description
%
%   [K,c] = conroaKstep(f,con,p,x,u,z,V,beta,gamma,s1,s2,sg)
%   [...] = conroaKstep(...,L1,L2,opts)
%
% Solves for a control function K which satisfies the SOS constraints 
% while holding all other variables fixed.
%
% Inputs:
%       -f:   polynomial vector field
%       -con: constraint function (scalar polynomial)
%       -p:   shape function (scalar polynomial)
%       -x:   state-space vector as PVAR
%       -u:   input vector as PVAR
%       -z:   Nz-by-1 column vector of monomials; specifies the control
%             function decision variable in the vector form K(x) = c'*z(x).
%       -V:   Lyapunov function (scalar polynomial)
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
%       -K: control function satisfying the SOS constraints; will be empty
%           if no feasible solution is found.
%       -c: Nz-by-1 vector of coefficients of the returned control
%           function.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2019-01-24
% * Changed:    2019-01-24
%
%% See also
%
% See PWROAEST, CONROAVSTEP
%%

if isempty(opts)
    opts = sosoptions;
end


% control decision variable
[K,c] = polydecvar('c',z);

% control system
fK = subs(f,u,K);

% control input constraint
cK = subs(con,u,K);

%% Constrained control feasibility problem
sosconstr = cell(2,1);

% {x: V(x) <= g} is contained in {x: grad(V)*fK < 0}
Vdot = jacobian(V,x)*fK;
sosconstr{1} = Vdot <= -L2 + s2*(V-gamma);

% {x: V(x) <= g} is contained in {x: cK(x) <= 0}
sosconstr{2} = -(cK + sg*(gamma-V)) >= 0;

% solve feasibility problem
sosconstr = vertcat(sosconstr{:});
[info,dopt] = sosopt(sosconstr,x,opts);

%% Output
if info.feas
    K = subs(K,dopt);
    c = subs(c,dopt);
else
    K=[];
    c=[];
end

end