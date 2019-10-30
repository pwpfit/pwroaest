function varargout = pwroaKstep(f1,f2,phi,con,~,x,u,z,V,~,gamma,~,s,si,sg,sj,~,L2,tau,roaopts)
% Solves the K-V-s step of the constraint-piecewise ROA iteration.
%
%% Usage & description
%
%   [K,c] = pwroavstep(f1,f2,phi,con,p,x,u,z,V,beta,gamma,s1,s2,si,sg,sj)
%   [...] = pwroavstep(...,L1,L2,opts)
%
% Solves for a common control function K which satisfies the SOS
% constraints while holding all other variables fixed.
%
% Inputs:
%       -f1:  first polynomial vector field
%       -f2:  second polynomial vector field
%       -phi: boundary condition (scalar field)
%       -con: constraint function (scalar polynomial)
%       -p:   shape function (scalar field)
%       -x:   state-space vector as PVAR
%       -u:   input vector as PVAR
%       -z:   Nz-by-1 column vector of monomials; specifies the Lyapunov
%             function decision variable in the vector form V(x) = c'*z(x).
%       -V:   Lyapunov function (scalar polynomial)
%       -beta:  level set of shape function
%       -gamma: level set of the Lyapunov function
%       -s1:  multiplier for beta-step; in the multiple V-step, 
%             multipliers [s1'; s2"].
%       -s2:  multipliers [s2'; s2"] for gamma-step (Lyapunov gradient)
%       -si:  multipliers [si'; si"] for gamma-step (boundary condition)
%       -sg:  multiplier for constraint step; in the multiple V-step,
%             multipliers [sg'; sg"] (Lyapunov gradient)
%       -sj:  multipliers [sj'; sj"] for constraint step (boundary
%             condition) in the multiple V-step.
%       -L1:  epsilon 1 (double or scalar field); enforces strict positive 
%             definiteness of the Lyapunov function.
%       -L2:  epsilon 2 (double or scalar field); enforces strict negative
%             definiteness of the Lyapunov gradients.
%       -tau: sample time
%       -opts:  SOS options structre; see SOSOPTIONS.
%
% Outputs:
%       -K: control function satisfying the SOS constraints; will be empty
%           if no feasible solution is found.
%       -c: Nz-by-1 vector of coefficients of the returned Lyapunov
%           function.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2018-05-23
% * Changed:    2019-10-28
%
%% See also
%
% See PWROAEST, ROAVSTEP
%%

if isa(roaopts,'roaoptions')
    opts = roaopts.sosopts;
else
    opts = roaopts;
end

if length(z) == 1
    % polynomial controller
    varargout = cell(1,2);
    
    % control decision variable
    [K,c] = polydecvar('c',z{1});

    % control system
    fK1 = subs(f1,u,K);
    fK2 = subs(f2,u,K);

    % control input constraint
    cK = subs(con,u,K);


    if length(V) == 1
        % common Lyapunov function
        V = V{1};

        %% Common K-V-s feasibility problem
        sosconstr = cell(3,1);

        % {x: V(x) <= g} is contained in {x: grad(V)*f < 0}
        [Vdot1,R1] = ddiff(V,x,tau,fK1);
        [Vdot2,R2] = ddiff(V,x,tau,fK2);

        % -( pa + (g-p2)*s - phi*si ) in SOS
        sosconstr{1} = -(Vdot1 + L2 + s(1)*(gamma-R1) - si(1)*phi) >= 0;
        % -( pb + (g-p2)*s + (phi-l)*si ) in SOS
        sosconstr{2} = -(Vdot2 + L2 + s(2)*(gamma-R2) + si(2)*(phi-L2)) >= 0;

        % {x: V(x) <= g} is contained in {x: c(x) <= 0}
        sosconstr{3} = -(cK + sg*(gamma-V)) >= 0;

        % solve problem
        sosconstr = vertcat(sosconstr{:});
        [info,dopt] = sosopt(sosconstr,x,opts);

    else
        % multiple Lyapunov functions
        V1 = V{1};
        V2 = V{2};

        %% Multiple K-V-s feasibility problem
        sosconstr = cell(4,1);

        % {x: V(x) <= g} is contained in {x: grad(V)*f < 0}
        [Vdot1,R1] = ddiff(V1,x,tau,fK1);
        [Vdot2,R2] = ddiff(V2,x,tau,fK2);
        % -( pa + (g-p2)*s - phi*si ) in SOS
        sosconstr{1} = -(Vdot1 + L2 + s(1)*(gamma-R1) - si(1)*phi) >= 0;
        % -( pb + (g-p2)*s + (phi-l)*si ) in SOS
        sosconstr{2} = -(Vdot2 + L2 + s(2)*(gamma-R2) + si(2)*(phi-L2)) >= 0;

        % {x: V1(x) <= g} intersects {x: phi(x) <= 0} is contained in {x: c(x) <= 0}
        sosconstr{3} = -(cK + sg(1)*(gamma - V1) - sj(:,1)*phi) >= 0;
        % {x: V2(x) <= g} intersects {x: phi(x) > 0} is contained in {x: c(x) <= 0}
        sosconstr{4} = -(cK + sg(2)*(gamma - V2) + sj(:,2)*(phi-L2)) >= 0;

        % solve problem
        sosconstr = vertcat(sosconstr{:});
        [info,dopt] = sosopt(sosconstr,x,opts);

    end

    %% Output
    if info.feas
        varargout{1} = subs(K,dopt);
        varargout{2} = subs(c,dopt);
    end
else
    % piecewise controller
    varargout = cell(1,4);

    % control decision variable
    [K1,c1] = polydecvar('c1',z{1});
    [K2,c2] = polydecvar('c2',z{2});

    % control system
    fK1 = subs(f1,u,K1);
    fK2 = subs(f2,u,K2);

    % control input constraint
    cK1 = subs(con,u,K1);
    cK2 = subs(con,u,K2);

    if length(V) == 1
        % common Lyapunov function
        V = V{1};

        %% Common K-V-s feasibility problem
        sosconstr = cell(3,1);

        % {x: V(x) <= g} is contained in {x: grad(V)*f < 0}
        [Vdot1,R1] = ddiff(V,x,tau,fK1);
        [Vdot2,R2] = ddiff(V,x,tau,fK2);

        % -( pa + (g-p2)*s - phi*si ) in SOS
        sosconstr{1} = -(Vdot1 + L2 + s(1)*(gamma-R1) - si(1)*phi) >= 0;
        % -( pb + (g-p2)*s + (phi-l)*si ) in SOS
        sosconstr{2} = -(Vdot2 + L2 + s(2)*(gamma-R2) + si(2)*(phi-L2)) >= 0;

        % {x: V1(x) <= g} intersects {x: phi(x) <= 0} is contained in {x: c(x) <= 0}
        sosconstr{3} = -(cK1 + sg(1)*(gamma - V) - sj(:,1)*phi) >= 0;
        % {x: V2(x) <= g} intersects {x: phi(x) > 0} is contained in {x: c(x) <= 0}
        sosconstr{4} = -(cK2 + sg(2)*(gamma - V) + sj(:,2)*(phi-L2)) >= 0;

        % solve problem
        sosconstr = vertcat(sosconstr{:});
        [info,dopt] = sosopt(sosconstr,x,opts);

    else
        % multiple Lyapunov functions
        V1 = V{1};
        V2 = V{2};

        %% Multiple K-V-s feasibility problem
        sosconstr = cell(4,1);

        % {x: V(x) <= g} is contained in {x: grad(V)*f < 0}
        [Vdot1,R1] = ddiff(V1,x,tau,fK1);
        [Vdot2,R2] = ddiff(V2,x,tau,fK2);
        % -( pa + (g-p2)*s - phi*si ) in SOS
        sosconstr{1} = -(Vdot1 + L2 + s(1)*(gamma-R1) - si(1)*phi) >= 0;
        % -( pb + (g-p2)*s + (phi-l)*si ) in SOS
        sosconstr{2} = -(Vdot2 + L2 + s(2)*(gamma-R2) + si(2)*(phi-L2)) >= 0;

        % {x: V1(x) <= g} intersects {x: phi(x) <= 0} is contained in {x: c(x) <= 0}
        sosconstr{3} = -(cK1 + sg(1)*(gamma - V1) - sj(:,1)*phi) >= 0;
        % {x: V2(x) <= g} intersects {x: phi(x) > 0} is contained in {x: c(x) <= 0}
        sosconstr{4} = -(cK2 + sg(2)*(gamma - V2) + sj(:,2)*(phi-L2)) >= 0;

        % solve problem
        sosconstr = vertcat(sosconstr{:});
        [info,dopt] = sosopt(sosconstr,x,opts);

    end

    %% Output
    if info.feas
        varargout{1} = subs(K1,dopt);
        varargout{2} = subs(K2,dopt);
        varargout{3} = subs(c1,dopt);
        varargout{4} = subs(c2,dopt);
    end
end