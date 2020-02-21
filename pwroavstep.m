function varargout = pwroavstep(f1,f2,phi,con,p,x,z,beta,gamma,s0,s,si,sg,sj,L1,L2,tau,roaopts)
% Solves the V-s step of the piecewise ROA iteration.
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
%       -con: constraint function (scalar polynomial)
%       -p:   shape function (scalar field)
%       -x:   state-space vector as PVAR
%       -z:   Nz-by-1 column vector of monomials; specifies the Lyapunov
%             function decision variable in the vector form V(x) = c'*z(x).
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
%        -tau: sample time
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
% * Changed:    2019-10-28
%
%% See also
%
% See PWROAEST, ROAVSTEP
%%

opts = roaopts.sosopts;
zi1  = roaopts.zi{1};
zi2  = roaopts.zi{end};
tau  = roaopts.tau;


if length(z) == 1
    % common Lyapunov function
    varargout = cell(1,2);
    
    % Lyapunov decision variable
    [V,c] = polydecvar('c',z{1});

    %% Common V-s feasibility problem
    sosconstr = cell(5,1);

    % V-L1 in SOS
    sosconstr{1} = V >= L1;

    % {x: p(x) <= b} is contained in {x: V(x) <= g}
    sosconstr{2} = -((V-gamma) + s0*(beta-p)) >= 0;

    % {x: V(x) <= g} is contained in {x: grad(V)*f < 0}
    [Vdot1,R1] = ddiff(V,x,tau,f1);
    [Vdot2,R2] = ddiff(V,x,tau,f2);
    % -( pa + (g-p2)*s - phi*si ) in SOS
    sosconstr{3} = -(Vdot1 + L2 + s(1)*(gamma-R1) - si(1)*phi) >= 0;
    % -( pb + (g-p2)*s + (phi-l)*si ) in SOS
    sosconstr{4} = -(Vdot2 + L2 + s(2)*(gamma-R2) + si(2)*(phi-L2)) >= 0;
    
    if length(con) == 1
        % polynomial control constraint
        con0 = con{1};
        
        % {x: V(x) <= g} is contained in {x: c(x) <= 0}
        sosconstr{5} = -(con0 + sg*(gamma-V)) >= 0;
    else
        % piecewise control constraint
        con1 = con{1};
        con2 = con{2};
        
        % {x: V(x) <= g} intersects {x: phi(x) <= 0} is contained in {x: c1(x) <= 0}
        sosconstr{7} = -(con1 + sg(1)*(gamma - V) - sj(:,1)*phi) >= 0;
        % {x: V(x) <= g} intersects {x: phi(x) > 0} is contained in {x: c2(x) <= 0}
        sosconstr{8} = -(con2 + sg(2)*(gamma - V) + sj(:,2)*(phi-L2)) >= 0;
    end

    % solve problem
    sosconstr = vertcat(sosconstr{:});
    [info,dopt] = sosopt(sosconstr,x,opts);

    % output
    if info.feas
        varargout{1}   = subs(V,dopt);
        varargout{end} = subs(c,dopt);
    end
else
    % multiple Lyapunov functions
    varargout = cell(1,4);

    % Lyapunov decision variables
    [V1,c1] = polydecvar('c1',z{1});
    [V2,c2] = polydecvar('c2',z{2});
    
    % control constraints
    con1 = con{1};
    con2 = con{end};
    
    % continuity decision variables
    ri = [ sosdecvar('ci1',zi1)
           sosdecvar('ci2',zi2)
           sosdecvar('ci3',zi1)
           sosdecvar('ci4',zi2)
    ];
    
    %% Multiple V-s feasibility problem
    sosconstr = cell(14,1);
    
    % Vi-L1 in SOS
    sosconstr{1} = V1 >= L1;
    sosconstr{2} = V2 >= L1;
    
    % {x: p(x) <= b} intersects {x: phi(x) <= 0} is contained in {x: V1(x) <= g}
    sosconstr{3} = -((V1-gamma) + s0(1)*(beta(1)-p(1))) >= 0;
    % {x: p(x) <= b} intersects {x: phi(x) > 0} is contained in {x: V2(x) <= g}
    sosconstr{4} = -((V2-gamma) + s0(2)*(beta(end)-p(end))) >= 0;
    
    % {x: V(x) <= g} is contained in {x: grad(V)*f < 0}
    [Vdot1,R1] = ddiff(V1,x,tau,f1);
    [Vdot2,R2] = ddiff(V2,x,tau,f2);
    % -( pa + (g-p2)*s - phi*si ) in SOS
    sosconstr{5} = -(Vdot1 + L2 + s(1)*(gamma-R1) - si(1)*phi) >= 0;
    % -( pb + (g-p2)*s + (phi-l)*si ) in SOS
    sosconstr{6} = -(Vdot2 + L2 + s(2)*(gamma-R2) + si(2)*(phi-L2)) >= 0;
    
    % {x: V1(x) <= g} intersects {x: phi(x) <= 0} is contained in {x: c(x) <= 0}
    sosconstr{7} = -(con1 + sg(1)*(gamma - V1) - sj(:,1)*phi) >= 0;
    % {x: V2(x) <= g} intersects {x: phi(x) > 0} is contained in {x: c(x) <= 0}
    sosconstr{8} = -(con2 + sg(2)*(gamma - V2) + sj(:,2)*(phi-L2)) >= 0;
    
    % {x: phi(x) <= 0} intersects {x: phi(x) >= 0} is contained in {x: V1(x) <= V2(x)}
    sosconstr{9} = -((V1-V2) + ri(1)*phi - ri(2)*phi) >= 0;
    % {x: phi(x) <= 0} intersects {x: phi(x) >= 0} is contained in {x: V1(x) >= V2(x)}
    sosconstr{10} = -((V2-V1) + ri(3)*phi - ri(4)*phi) >= 0;
    
    sosconstr{11} = ri(1) >= 0;
    sosconstr{12} = ri(2) >= 0;
    sosconstr{13} = ri(3) >= 0;
    sosconstr{14} = ri(4) >= 0;
    
    % solve problem
    sosconstr = vertcat(sosconstr{:});
    [info,dopt] = sosopt(sosconstr,x,opts);

    % output
    if info.feas
        varargout{1} = subs(V1,dopt);
        varargout{2} = subs(V2,dopt);
        varargout{3} = subs(c1,dopt);
        varargout{4} = subs(c2,dopt);
    end

end    