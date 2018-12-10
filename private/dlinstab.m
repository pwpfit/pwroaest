function [V,A,P]=dlinstab(f,x,Q,tau)
% Linear stability analysis for discrete system, x+ = x + tau*f
if nargin==2 || isempty(Q)
    Nx = length(x);
    Q=eye(Nx);
end

if tau == 0
    % fall back to continuous analysis
    [V,A,P] = linstab(f,x,Q);
    return
elseif tau > 0
    % Linearize: x+ = A*x
    A=plinearize(x+tau*f,x);
else
    % implicit discretization
    % Linearize: x- = A^-1*x
    A=plinearize(x+tau*f,x)^-1;
end        

% If A is stable then solve the Lyapunov equation: A'*P+P*A = -Q
ev = eig(A);
if max(abs(ev)) < 1
    I = eye(size(A));
    if exist('dlyap','file')==2
        P=dlyap(A',Q);
    else
        %P=LOCALlyap(A',Q);
        error('Discrete Lyapunov equation not found.');
    end
    V=x'*P*x;
else
    P=[];
    V=[];
end

