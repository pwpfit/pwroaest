function [A,B,f0] = dplinearize(f,x,u,tau,x0,u0)
% Linearize polynomial function f+(x,u) = x + tau*f(x,u) around trim
% condition x = x0 and u = u0.

% Error Checking
if nargin==2
    %   [A,f0] = dplinearize(f,x)
    x0 = [];
    u0 = [];
    tau = 0;
    u = [];
elseif nargin==3
    if ispvar(u)        
        %   [A,B,f0] = dplinearize(f,x,u)
        x0 = [];
        u0 = [];
        tau = 0;
    else
        %   [A,f0] = dplinearize(f,x,tau)
        x0 = [];
        u0 = [];
        tau = u;
        u = [];
    end
elseif nargin==4
    if ispvar(u)
        %   [A,B,f0] = dplinearize(f,x,u,tau)
        x0 = [];
        u0 = [];
    elseif isempty(u)
        %   [A,f0] = dplinearize(f,x,[],x0)
        x0 = tau;
        u0 = [];
        tau = 0;
        u = [];
    else
        %   [A,f0] = dplinearize(f,x,tau,x0)
        x0 = tau;
        u0 = [];
        tau = u;
        u = [];
    end    
elseif nargin==5
    if isempty(tau)
        %   [A,B,f0] = dplinearize(f,x,u,[],x0)
        tau = 0;
    end
    % else: [A,B,f0] = dplinearize(f,x,u,tau,x0)
    u0 = [];
elseif isempty(tau)
    %   [A,B,f0] = dplinearize(f,x,u,[],x0,u0)
    tau = 0;
end

if tau == 0
    % fall back to continuous linearization
    % xdot = f0 + A*(x-x0) + B*(u-u0) + O(x,u)
    [A,B,f0] = plinearize(f,x,u,x0,u0);
    return
end

% else:
% linearize x+ = f0 + A*(x-x0) + B*(u-u0) + O(x,u)
[A,B,f0] = plinearize(x+tau*f,x,u,x0,u0);

% Implicit discretization
if tau < 0
    % x- = -A^-1*f0 + A^-1*x - A^-1*B*u - A^-1*O(x,u)
    f0 = -A\f0;
    B = -A\B;
    A = A^-1;
end
