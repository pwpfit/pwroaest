function [beta,V,gamma,K,iter] = pwroaest(f1,f2,phi,x,roaopts)
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
% * Changed:    2019-01-15
%
%% See also
%
% See ROAEST, PWROAOPTIONS, PCONTAIN, PWPCONTAIN
%%

% information from options
p  = roaopts.p;
zV = roaopts.zVi;
z1 = roaopts.z1i;
z2 = roaopts.z2i;
zg = roaopts.zgi;
zi = roaopts.zi;
zK = roaopts.zKi;
L2 = roaopts.L2;
L1 = roaopts.L1;
Q  = roaopts.Q;
c  = roaopts.c;

NstepBis = roaopts.NstepBis;
sopts = roaopts.sosopts;
gopts = roaopts.gsosopts;
gopts.minobj = 0;
gammamax = roaopts.gammamax;
betamax = roaopts.betamax;
Vin = roaopts.Vi0;
Kin = roaopts.Ki0;

u  = roaopts.u;

display = roaopts.display;
debug   = roaopts.debug;
log = roaopts.log;
logpath = roaopts.logpath;

if ~strcmp(log,'none') && ~exist(logpath{1},'file')
    [success,info] = mkdir(logpath{1});
    if success && ~isempty(info)
        warning(info)
    elseif ~success
        error(info)
    end
elseif exist(logpath{1},'file')
    delete([logpath{1} '/*.mat']);
end

% Vdeg = zV.maxdeg;
Nsteps = NstepBis;

% initialize storage
c0 = cell(Nsteps,1);
iter= struct('V',c0,'K',c0,'beta',c0,'gamma',c0,'s0',c0,'s',c0,'si',c0,'sg',c0,'sj',c0,'time',c0,'aux',c0);

% boundary condition at origin
phi0 = double(subs(phi, x, zeros(size(x))));

% Lyapunov functions
V = cell(size(zV));

% controllers
K  = cell(size(zK));
cK = cell(size(zK));

%% Run V-s iteration
fprintf('\n---------------Beginning piecewise V-s iteration\n');
biscount = 0;
for i1=1:NstepBis
    biscount = biscount+1;

    tic;

    if i1==1
        if ~isempty(Kin)
            K(:) = Kin;
        else
            K(:) = zeros(size(u));
        end
    elseif isempty(zK)
        % nothing to do
        
    elseif gpre <= gmin
        % local K-V-s problem
        K{1} = roaKstep(f1,c,p,x,u,zK{1},V{1},b,g,s0,s,sg,L1,L2,sopts);
        if isempty(K{1})
            if strcmp(display,'on')
                fprintf('local K-step infeasible at iteration = %d\n',i1);
            end
            break;
        elseif length(K) > 1
            K{2} = K{1};
        end
    else
        [K{:}] = pwroaKstep(f1,f2,phi,c,p,x,u,zK,V,[b1 b2],g,s0,s,si,sg,sj,L1,L2,roaopts);
        if isempty(K{1})
            if strcmp(display,'on')
                fprintf('common K-step infeasible at iteration = %d\n',i1);
            end
            break;
        end
    end        
    
    % system & constraints under control
    if ~isempty(u)
        fK1 = subs(f1,u,K{1});
        fK2 = subs(f2,u,K{end});
        cK{1}   = polynomial(subs(c,u,K{1}));
        cK{end} = polynomial(subs(c,u,K{end}));
    else
        fK1 = f1;
        fK2 = f2;
        cK  = {c};
    end
    
    %======================================================================
    % Find V step:
    % Hold s1, s2, si, b, g fixed and solve the Gamma Step and Beta Step
    % for Appropiate Lyapunov function with an additional constraint
    % V-L1 in SOS
    %======================================================================
    if i1==1
        if ~isempty(Vin)
            V = Vin;
        elseif phi0 < 0
            % Construct Lyap function from linearization of LHS function
            V{1}=linstab(fK1,x,Q);
            if length(V) > 1
                V{2} = V{1};
            end
        else
            [V{:}]=pwlinstab(fK1,fK2,phi,x);
        end
        
    elseif gpre <= gmin
        % local V-s problem
        [V{1},~] = conroavstep(fK1,cK{1},p,x,zV{1},b,g,s0,s,sg,L1,L2,sopts);
        if isempty(V{1})
            if strcmp(display,'on')
                fprintf('local V-step infeasible at iteration = %d\n',i1);
            end
            break;
        elseif length(V) > 1
            V{2} = V{1};
        end
    else
        [V{:}] = pwroavstep(fK1,fK2,phi,cK,p,x,zV,[b1 b2],g,s0,s,si,sg,sj,L1,L2,roaopts);
        if isempty(V{1})
            if strcmp(display,'on') && length(V) == 1
                fprintf('common V-step infeasible at iteration = %d\n',i1);
            elseif strcmp(display,'on')
                fprintf('multiple V-step infeasible at iteration = %d\n',i1);
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
        [gbnds,s] = pcontain(jacobian(V{1},x)*fK1+L2,V{1},z2{1},gopts);
        if isempty(gbnds)
            if strcmp(display,'on')
                fprintf('pre gamma step infeasible at iteration = %d\n',i1);
            end
            break;
        end
        gpre = gbnds(1);
        
        %======================================================================
        % Min Gamma Step: Solve the problem
        % {x:V(x) <= gamma} is contained in {x:phi(x)<=0}
        %======================================================================
        gopts.maxobj = gammamax;
        [gbnds,~] = pcontain(phi,V{1},[],gopts);
        if isempty(gbnds)
            if strcmp(display,'on')
                fprintf('min gamma step infeasible at iteration = %d\n',i1);
            end
            break;
        end
        gmin = gbnds(1);
        
        if strcmp(debug,'on')
            fprintf('debug: gpre = %4.6f \t gmin = %4.6f\n', gpre, gmin);
        end

    else
        % origin at boundary
        % no local problem possible
        gpre = [];
        gmin = [];
    end
    
    if gpre <= gmin
        % estimated region of attraction does not reach boundary
        gstb = gpre;
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
        [gbnds,s1,si1] = pwpcontain(jacobian(V{1},x)*fK1+L2,  ...
                                    V{1}, phi, z2{1}, zi{1}, gopts ...
        );
        if isempty(gbnds)
            if strcmp(display,'on')
                fprintf('gamma 1 step infeasible at iteration = %d\n',i1);
            end
            break;
        end
        g1 = gbnds(1);

        %==================================================================
        % Gamma 2 Step: Solve the following problems
        % {x:V(x) <= gamma2} intersects {x:phi(x) > 0} is contained 
        % in {x:grad(V)*f2 < 0}
        %
        % max gamma2 subject to
        %     -[grad(V)*f2 + (gamma - V)*s + phi*si] in SOS,  s,si in SOS
        %==================================================================
        gopts.maxobj = gammamax;
        [gbnds,s2,si2] = pwpcontain(jacobian(V{end},x)*fK2+L2,  ...
                                    V{end}, -phi+L2, z2{end}, zi{end}, gopts ...
        );
        if isempty(gbnds)
            if strcmp(display,'on')
                fprintf('gamma 2 step infeasible at iteration = %d\n',i1);
            end
            break;
        end
        g2 = gbnds(1);

        if strcmp(debug,'on')
            fprintf('debug: g1 = %4.6f \t g2 = %4.6f\n', g1, g2);
        end
        
        gstb  = min(g1,g2);
    
        s  = [s1  s2 ];
        si = [si1 si2];
    end

    if gstb > .99*gammamax
        if strcmp(display,'on')
            fprintf('result of gamma step close to maximum (99%%) at iteration = %d\n',i1);
        end
        break;
    end 
    
    if (length(V) == 1 && length(cK) == 1) || any(gpre <= gmin)
        %==================================================================
        % Constraint Step: Solve the following problem
        % {x: V(x) <= gamma} is contained in {x: c(x) <= 0}
        %==================================================================
        [gcons,sg] = pcontain(cK{1},V{1},zg{1},gopts);
        if isempty(gcons)
            if strcmp(display,'on')
                fprintf('constraint step infeasible at iteration = %d\n',i1);
            end
            break;
        end

        gcon = gcons(1);
        sj = polynomial;
        
        if strcmp(debug,'on')
            fprintf('debug: gstb = %4.6f \t gcon = %4.6f\n', gstb, gcon);
        end
        
        g = min([gstb, gcon]);
        
        %======================================================================
        % Beta Step: Solve the following problem
        % {x: p(x)) <= beta} is contained in {x: V(x) <= gamma}
        % max beta subject to
        %                 -[(V - gamma) + (beta - p)*s0] in SOS, s0 in SOS
        %======================================================================
        gopts.maxobj = betamax;
        [bbnds,s0]=pcontain(V{1}-g,p,z1{1},gopts);
        if isempty(bbnds)
            if strcmp(display,'on')
                fprintf('beta step infeasible at iteration = %d\n',i1);
            end
            break;
        end
        b1 = bbnds(1);
        
        b2 = [];
    else
        %==================================================================
        % Constraint 1 Step: Solve the following problem
        % {x: V1(x) <= gamma} intersects {x:phi(x) <= 0} is contained in 
        % {x: c1(x) <= 0}
        %==================================================================
        [gcons,sg1,sj1] = pwpcontain(cK{1},V{1},phi,zg{1},zi{1},gopts);
        if isempty(gcons)
            if strcmp(display,'on')
                fprintf('constraint 1 step infeasible at iteration = %d\n',i1);
            end
            break;
        end

        gc1 = gcons(1);
        
        %==================================================================
        % Constraint 2 Step: Solve the following problem
        % {x: V2(x) <= gamma} intersects {x:phi(x) > 0} is contained in 
        % {x: c2(x) <= 0}
        %==================================================================
        [gcons,sg2,sj2] = pwpcontain(cK{end},V{end},-phi+L2,zg{end},zi{end},gopts);
        if isempty(gcons)
            if strcmp(display,'on')
                fprintf('constraint 2 step infeasible at iteration = %d\n',i1);
            end
            break;
        end
        
        gc2 = gcons(1);
        
        if strcmp(debug,'on')
            fprintf('debug: gstb = %4.6f \t gc1 = %4.6f \t gc2 = %4.6f\n', gstb, gc1, gc2);
        end
        
        gcon = min(gc1,gc2);
        g = min([gstb, gcon]);
        
        sg = [sg1 sg2];
        sj = [sj1 sj2];
        
        %======================================================================
        % Beta 1 Step: Solve the following problem
        % {x: p(x)) <= beta1} intersects {x:phi(x) <= 0} is contained 
        % in {x: V(x) <= gamma}
        % max beta subject to
        %   -[(V - gamma1) + (beta1 - p)*s0 - phi*si] in SOS, s0,si in SOS
        %======================================================================
        gopts.maxobj = betamax;
%         [bbnds,sb1,sj1]=pwpcontain(V{1}-g1, ...
%                                    p, phi, z1, zi, gopts ...
%         );
        [bbnds,sb1]=pcontain(V{1}-g,p,z1{1},gopts);
        if isempty(bbnds)
            if strcmp(display,'on')
                fprintf('beta 1 step infeasible at iteration = %d\n',i1);
            end
            break;
        end
        b1 = bbnds(1);
        
        %======================================================================
        % Beta 2 Step: Solve the following problem
        % {x: p(x)) <= beta2} intersects {x:phi(x) > 0} is contained 
        % in {x: V(x) <= gamma}
        % max beta2 subject to
        %   -[(V - gamma2) + (beta2 - p)*s0 + phi*si] in SOS, s0,si in SOS
        %======================================================================
        gopts.maxobj = betamax;
%         [bbnds,sb2,sj2]=pwpcontain(V{2}-g2, ...
%                                    p, -phi+L2, z1, zi, gopts ...
%         );
        [bbnds,sb2]=pcontain(V{end}-g,p,z1{end},gopts);
        if isempty(bbnds)
            if strcmp(display,'on')
                fprintf('beta 2 step infeasible at iteration = %d\n',i1);
            end
            break;
        end
        b2 = bbnds(1);
        
        if strcmp(debug,'on')
            fprintf('debug: b1 = %4.6f \t b2 = %4.6f\n', b1, b2);
        end
        
        s0 = [sb1 sb2];
%         sj = [0 0]; %sj = [sj1 sj2];
    end

    b  = min([b1 b2]);

    if b > .99*betamax && strcmp(debug,'on')
        fprintf('warning: result of beta step close to maximum (99%%) at iteration = %d\n',i1);
    end
    
    
    if ~strcmp(roaopts.gammacheck, 'none') && gpre > gmin
        %==============================================================
        % Gamma Check Step: Solve the following problems
        % {x:V(x) <= gamma} intersects {x:phi(x) <= 0} is contained 
        % in {x:grad(V)*f1 < 0}
        %
        % and
        %
        % {x:V(x) <= gamma} intersects {x:phi(x) > 0} is contained 
        % in {x:grad(V)*f2 < 0}
        %
        % with gamma = min(gamma1,gamma2)
        %==============================================================
        [s1,si1] = pwpcontain_check(jacobian(V{1},x)*fK1+L2,  ...
                                    V{1}-g, phi, z2{1}, zi{1}, gopts ...
        );
        [s2,si2] = pwpcontain_check(jacobian(V{end},x)*fK2+L2,  ...
                                    V{end}-g, -phi+L2, z2{end}, zi{end}, gopts ...
        );

        if isempty(s1) || isempty(s2)
            if strcmp(display,'on')
                fprintf('gamma check step infeasible at iteration = %d\n',i1);
            end
            if strcmp(roaopts.gammacheck, 'check')
                break;
            end
        else    
            s  = [s1  s2 ];
            si = [si1 si2];
        end
    end
    
    % Print results and store iteration data
    if strcmp(display,'on')
        fprintf('iteration = %d  \t beta = %4.6f \t gamma = %4.6f\n',i1,b,g);
    end
    iteration.V     = V;
    iteration.K     = K;
    iteration.beta  = b;
    iteration.gamma = g;
    iteration.s0    = s0;
    iteration.s     = s;
    iteration.si    = si;
    iteration.sg    = sg;
    iteration.sj    = sj;
    iteration.aux   = struct('gstb',gstb,'gcon',gcon,'gpre',gpre,'gmin',gmin);
    iteration.time  = toc;
    iter(i1) = iteration;
    if strcmp(log,'step')
        save([logpath{:} sprintf('iter%d',i1)], '-struct', 'iteration');
    end
end
if strcmp(display,'on')
    fprintf('---------------Ending V-s iteration.\n');
end
    
%% Outputs
iter(biscount+1:end) = [];
[~, idx] = max([iter.beta -1]);     % handle empty beta value(s)
result = iter(idx);
beta  = result.beta;
V     = result.V;
K     = result.K;
gamma = result.gamma;   % handle empty gamma value


if any(strcmp(log,{'step' 'result'}))
    save([logpath{:} 'result'], 'iter');
    save([logpath{:} 'result'], '-struct', 'result', '-append');
end
