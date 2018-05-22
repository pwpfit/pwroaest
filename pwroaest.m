function [beta,V,gamma,s1,s2,si,iter] = pwroaest(f1,f2,x,phi,roaopts)
% Estimates lower bound of piece-wise region of attraction.
%
%% Usage & description
%
%   [beta,V,gamma,s1,s2,si,iter] = pwroaest(f1, f2, x, phi, ropts)
%
% Estimates the lower bound of the region-of-attraction of the piecewise
% polynomial vector field
%
%          / f1(x)      if phi(x) <= 0;
%   xdot = |
%          \ f2(x)      else;
%
% about the equilibrium point x = 0 with phi(0) < 0.
%
% The piecewise ROA estimation problem is formulated as:
%   max beta
%   subject to:
%         V-L1 in SOS
%       -( (V-gamma) + (beta-p)*s1 ) in SOS
%       -( (grad(V)*f1 + L2) + (gamma-V)*s2' - phi*si' ) in SOS
%       -( (grad(V)*f2 + L2) + (gamma-V)*s2" + phi*si" ) in SOS
%
% The problem above is solved by a V-s iteration algorithm based on the
% work of Balas et al. (2009).
%
% Inputs:
%       -f1:  first polynomial vector field
%       -f2:  second polynomial vector field
%       -x:   state-space vector as PVAR
%       -phi: boundary condition (scalar field)
%       -ropts: options for piece-wise ROA estimation; see PWROAOPTIONS.
%
% Outputs:
%       -beta: maximum beta
%       -V:  Lyapunov function corresponding to beta
%       -s1: multiplier for beta-step 
%       -s2: multipliers [s2'; s2"] for gamma-step (Lyapunov gradient)
%       -si: multipliers [si'; si"] for gamma-step (boundary condition)
%       -iter: structure with fields V, beta, gamma, s1, s2, and time;
%              contains the results for each iteration.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2018-05-22
% * Changed:    2018-05-22
%
%% See also
%
% See ROAEST, PWROAOPTIONS
%%

% information from options
p  = roaopts.p;
zV = roaopts.zV;
z1 = roaopts.z1;
z2 = roaopts.z2;
zi = roaopts.zi;
L2 = roaopts.L2;
L1 = roaopts.L1;
Q  = roaopts.Q;
NstepBis = roaopts.NstepBis;
sopts = roaopts.sosopts;
gopts = roaopts.gsosopts;
gopts.minobj = 0;
gammamax = roaopts.gammamax;
betamax = roaopts.betamax;
display = roaopts.display;
Vin = roaopts.Vin;

Vdeg = zV.maxdeg;
Nsteps = NstepBis;

% initialize storage
c0 = cell(Nsteps,1);
iter= struct('V',c0,'beta',c0,'gamma',c0,'s1',c0,'s2',c0,'time',c0);


%% Run V-s iteration
fprintf('\n---------------Beginning piecewise V-s iteration\n');
biscount = 0;
for i1=1:NstepBis
    tic;
    
    %======================================================================
    % Find V step:
    % Hold s1, s2, si, b, g fixed and solve the Gamma Step and Beta Step
    % for Appropiate Lyapunov function with an additional constraint
    % V-L1 in SOS
    %======================================================================
    if i1==1
        if isempty(Vin)
            % Construct Lyap function from linearization
            V=linstab(f1,x,Q);
        else
            V = Vin;
        end
    else
%         [V,~] = roavstep(f,p,x,zV,b,g,s1,s2,L1,L2,sopts);
        V = [];
        if isempty(V)
            if strcmp(display,'on')
                fprintf('V-step infeasible at iteration = %d\n',i1);
            end
            break;
        end
    end

    %======================================================================
    % Gamma Step: Solve the following problem
    % {x:V(x) <= gamma} intersects {x:phi(x) <= 0} is contained 
    % in {x:grad(V)*f1 < 0}
    % AND
    % {x:V(x) <= gamma} intersects {x:phi(x) >= 0} is contained 
    % in {x:grad(V)*f2 < 0}
    %
    % max gamma subject to
    %     -[grad(V)*f1/2 + (gamma - V)*s2 -/+ phi*si] in SOS, s2, si in SOS
    %======================================================================
    gopts.maxobj = gammamax;
%     [gbnds,s2]=pcontain(jacobian(V,x)*f+L2,V,z2,gopts);
    [gbnds,s2,si] = pwpcontain(jacobian(V,x)*f1+L2,  ...
                               jacobian(V,x)*f2+L2,  ...
                               V, phi, z2, zi, L2, gopts ...
	);
    if isempty(gbnds)
        if strcmp(display,'on')
            fprintf('gamma step infeasible at iteration = %d\n',i1);
        end
        break;
    end
    g = gbnds(1);

    %======================================================================
    % Beta Step: Solve the following problem
    % {x: p(x)) <= beta} is contained in {x: V(x) <= gamma}
    % max beta subject to
    %                 -[(V - gamma) + (beta - p)*s1] in SOS, s1 in SOS
    %======================================================================
    gopts.maxobj = betamax;
    [bbnds,s1]=pcontain(V-g,p,z1,gopts);
    if isempty(bbnds)
        if strcmp(display,'on')
            fprintf('beta step infeasible at iteration = %d\n',i1);
        end
        break;
    end
    b = bbnds(1);
    
    % Print results and store iteration data
    if strcmp(display,'on')
        fprintf('iteration = %d  \t beta = %4.6f \t gamma = %4.6f\n',i1,b,g);
    end
    iter(i1).V     = V;
    iter(i1).beta  = b;
    iter(i1).gamma = g;
    iter(i1).s1    = s1;
    iter(i1).s2    = s2;
    iter(i1).si    = si;
    iter(i1).time  = toc;
    biscount = biscount+1;
end
if strcmp(display,'on')
    fprintf('---------------Ending V-s iteration.\n');
end
    
%% Outputs
[~, idx] = max([iter.beta]);
beta  = iter(idx).beta;
V     = iter(idx).V;
s1    = iter(idx).s1;
s2    = iter(idx).s2;
gamma = iter(idx).gamma;

end