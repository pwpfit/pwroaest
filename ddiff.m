function [fdot,f] = ddiff(f,v,tau,vdot)
% Discrete derivative, ddiff(f) = f(v+) - f(v).

if tau == 0
    % fall back to continuous derivative
    fdot = jacobian(f,v)*vdot;
elseif tau > 0
    % explicit discretization
    % v+ = v + tau*vdot
    fdot = subs(f,v,v+tau*vdot) - f;
else
    % implicit discretization
    % v- = v - |tau|*vdot
    fdot = f - subs(f,v,v+tau*vdot);
    f    = subs(f,v,v+tau*vdot);
end

end