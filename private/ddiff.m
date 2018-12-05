function fdot = ddiff(f,v,tau,vdot)
% Discrete derivative, ddiff(f) = f(v+tau*vdot) - f(v).

if tau == 0
    % fall back to continuous derivative
    fdot = jacobian(f,v)*vdot;
else
    fdot = subs(f,v,v+tau*vdot)/tau - f;
end

end