function [beta,V,gamma,varargout] = pwroaest(f1,f2,phi,x,roaopts)
% Estimates lower bound of piece-wise region of attraction.
%
%% Usage & description
%
%   [beta,V,gamma,iter] = pwroaest(f1, f2, phi, x, ropts)
%   [beta,V,gamma,s1,s2,si,iter] = pwroaest(...)
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
%       -phi: boundary condition (scalar field)
%       -x:   state-space vector as PVAR
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
% * Changed:    2018-05-23
%
%% See also
%
% See ROAEST, PWROAOPTIONS, PCONTAIN, PWPCONTAIN
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
iter= struct('V',c0,'beta',c0,'gamma',c0,'s0',c0,'s',c0,'si',c0,'time',c0);

% boundary condition at origin
phi0 = double(subs(phi, x, zeros(size(x))));


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
        g = 1;
        
        if ~isempty(Vin)
            V = Vin;
        elseif phi0 < 0
            % Construct Lyap function from linearization of LHS function
            V=linstab(f1,x,Q);
        else
            V=pwlinstab(f1,f2,phi,x);
        end
    elseif g <= gmin
        % local V-s problem
        [V,~] = roavstep(f1,p,x,zV,b,g,s0,s,L1,L2,sopts);
        if isempty(V)
            if strcmp(display,'on')
                fprintf('local V-step infeasible at iteration = %d\n',i1);
            end
            break;
        end
    else
        [V,~] = pwroavstep(f1,f2,phi,p,x,zV,b,g,s0,s,si,L1,L2,sopts);
        if isempty(V)
            if strcmp(display,'on')
                fprintf('common V-step infeasible at iteration = %d\n',i1);
            end
            break;
        end
    end
    
    if phi0 < 0
    %======================================================================
    % Pre Gamma Step: Solve the problem
    % {x:V(x) <= gamma} is contained in {x:grad(V)*f1 < 0}
    %======================================================================
    gopts.maxobj = gammamax;
    [gbnds,s] = pcontain(jacobian(V,x)*f1+L2,V,z2,gopts);
    if isempty(gbnds)
        if strcmp(display,'on')
            fprintf('pre gamma step infeasible at iteration = %d\n',i1);
        end
        break;
    end
    gpre = gbnds(1)
    
    %======================================================================
    % Min Gamma Step: Solve the problem
    % {x:V(x) <= gamma} is contained in {x:phi(x)<=0}
    %======================================================================
    gopts.maxobj = gammamax;
    [gbnds,~] = pcontain(phi,V,[],gopts);
    if isempty(gbnds)
        if strcmp(display,'on')
            fprintf('min gamma step infeasible at iteration = %d\n',i1);
        end
        break;
    end
    gmin = gbnds(1)
    
    else
        % origin at boundary
        % no local problem possible
        gpre = [];
        gmin = [];
    end
    
    if gpre <= gmin
        % estimated region of attraction does not reach boundary
        g  = gpre;
        si = polynomial;
        if strcmp(display,'on')
            fprintf('Local polynomial problem at iteration = %d\n',i1);
        end
    else
        %==================================================================
        % Gamma 1 Step: Solve the following problems
        % {x:V(x) <= gamma1} intersects {x:phi(x) <= 0} is contained 
        % in {x:grad(V)*f1 < 0}
        % 
        % max gamma1 subject to
        %     -[grad(V)*f1 + (gamma1 - V)*s - phi*si] in SOS,  s,si in SOS
        %==================================================================
        gopts.maxobj = gammamax;
        [gbnds,s1,si1] = pwpcontain(jacobian(V,x)*f1+L2,  ...
                                    V, phi, z2, zi, gopts ...
        );
        if isempty(gbnds)
            if strcmp(display,'on')
                fprintf('gamma 1 step infeasible at iteration = %d\n',i1);
            end
            break;
        end
        g1 = gbnds(1)

        %==================================================================
        % Gamma 2 Step: Solve the following problems
        % {x:V(x) <= gamma2} intersects {x:phi(x) > 0} is contained 
        % in {x:grad(V)*f2 < 0}
        %
        % max gamma2 subject to
        %     -[grad(V)*f2 + (gamma - V)*s + phi*si] in SOS,  s,si in SOS
        %==================================================================
        gopts.maxobj = gammamax;
        [gbnds,s2,si2] = pwpcontain(jacobian(V,x)*f2+L2,  ...
                                    V, -phi+L2, z2, zi, gopts ...
        );
        if isempty(gbnds)
            if strcmp(display,'on')
                fprintf('gamma 2 step infeasible at iteration = %d\n',i1);
            end
            break;
        end
        g2 = gbnds(1)

        s  = [s1  s2 ];
        si = [si1 si2];
        g  = min(g1,g2);
    end

    if g > .99*gammamax
        if strcmp(display,'on')
            fprintf('result of gamma step close to maximum (99%%) at iteration = %d\n',i1);
        end
        break;
    end 
    
    %======================================================================
    % Beta Step: Solve the following problem
    % {x: p(x)) <= beta} is contained in {x: V(x) <= gamma}
    % max beta subject to
    %                 -[(V - gamma) + (beta - p)*s1] in SOS, s1 in SOS
    %======================================================================
    gopts.maxobj = betamax;
    [bbnds,s0]=pcontain(V-g,p,z1,gopts);
    if isempty(bbnds)
        if strcmp(display,'on')
            fprintf('beta step infeasible at iteration = %d\n',i1);
        end
        break;
    end
    b = bbnds(1);
    
    if b > .99*betamax && strcmp(display,'on')
        fprintf('warning: result of beta step close to maximum (99%%) at iteration = %d\n',i1);
    end
    
    % Print results and store iteration data
    if strcmp(display,'on')
        fprintf('iteration = %d  \t beta = %4.6f \t gamma = %4.6f\n',i1,b,g);
    end
    iter(i1).V     = V;
    iter(i1).beta  = b;
    iter(i1).gamma = [g gpre gmin];
    iter(i1).s0    = s0;
    iter(i1).s     = s;
    iter(i1).si    = si;
    iter(i1).time  = toc;
    biscount = biscount+1;
end
if strcmp(display,'on')
    fprintf('---------------Ending V-s iteration.\n');
end
    
%% Outputs
iter(biscount+1:end) = [];
[~, idx] = max([iter.beta]);
beta  = iter(idx).beta;
V     = iter(idx).V;
s0    = iter(idx).s0;
s2    = iter(idx).s;
si    = iter(idx).si;
gamma = iter(idx).gamma(1);

if nargout <= 4
    varargout = {iter};
else
    varargout = {s0,s2,si,iter};

end