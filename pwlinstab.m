function [V, A1, A2, P] = pwlinstab(f1, f2, phi, x, Q)

if ~exist('Q', 'var')
    Q = 1e-6;
end

% Linearize: xdot = A*x
A1 = plinearize(f1,x);
A2 = plinearize(f2,x);

ev1 = eig(A1);
ev2 = eig(A2);
if min(real(ev1), real(ev2)) < 0
    % if A1, A2 are stable solve LMI
    %   A1'*P + P*A1 = -Q
    %   A2'*P + P*A2 = -Q
    [V,P] = sosdecvar('p', x);
    
    sosc = polyconstr;
    sosc(1) = x'*(A1'*P + P*A1)*x - phi <= -x'*Q*x;
    sosc(2) = x'*(A2'*P + P*A2)*x + phi <= -x'*Q*x;
    sosc(3) = V >= x'*Q*x;
    
    [info,dopt] = sosopt(sosc, x, sosoptions);
else
    info.feas = false;
end

if info.feas
    P = subs(P, dopt);
    V = x'*P*x;
else
    V = [];
    P = [];
end