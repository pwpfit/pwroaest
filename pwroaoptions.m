classdef pwroaoptions < conroaoptions
% Options for piecewise ROA estimation.
%
%% Usage & description
%
%   opts = pwroaoptions(f1, f2, phi, x)
%
% Standard options for piecewise ROA estimation;
%
%   opts = pwroaoptions(f, x)
%
% Fall-back to single polynomial ROA estimation;
%
%   opts = pwroaoptions(..., u)
%
% Optional specification of control input;
%
%   opts = pwroaoptions(..., 'Name', Value, ...)
%
% Specifiy list of name-value pairs for additional options (see below).
%
% with inputs
%       -f,f1:  first (single) polynomial vector field
%       -f2:  second polynomial vector field
%       -phi: boundary condition (scalar field)
%       -x:   state-space vector as PVAR
%       -u:   optional input vector as PVAR
%
%   Name, Value
%       Basic options:
%       -c:     constraints as polynomial vector field
%       -Kin:   Initial control function [empty by default]
%       -Ki0:   Cell of initial control functions [empty by default]
%       -L1:    Margin to enforce positive definiteness of the Lyapunov
%               function [default = 0]
%       -L2:    Margin to enforce positive definiteness of the gradient of
%               the Lyapunov function [default = 1e-6*x^2]
%       -p:     Ellipsoidal shape [default = x'*x]
%               If p is an empty polynomial, a set-inclusion constraint is
%               implemented in the V-step; here, the next Lyapunov function
%               is constrained to include the previous region-of-attraction
%               estimated rather than an ellipsoid shape. For comparison,
%               the beta step is executed and logged for the unit ball but 
%               not counted towards execution time.
%       -Q:     Positive definite matrix for initial linearized Lyapunov
%               equation [see pwlinstab, default = 1]
%       -tau:   sample time for discrete systems [default = 0]
%       -Vin:   Initial Lyapunov function [empty by default]
%       -Vi0:   Cell of initial multiple Lyapunov functions [empty by
%               default]
%       -u:  Control input vector (see above)
%       -xi: Variables of boundary condition as PVAR; must be subset of the
%            state-space variables. [default = opts.x] DEPRECATED
%
%       Monomial vectors:
%       -z1: Monomials for shape multiplier in beta-s1 step
%       -z2: Monomials for gradient multiplier in gamma-s2-[si] step
%       -zi: Monomials for boundary multiplier in gamma-s2-si step
%       -zg: Monomials for constraint multiplier
%            Specifically, z* are column vectors of monomials used to 
%            specify s*(x) in the Gram matrix form, s*(x)=z*(x)'*C*z*(x). 
%            [default =  monomials(x, 1:2)]
%
%       -zK: Vector of monomials for control function;
%       -zV: Vector of monomials for Lyapunov function;
%            Specifically, zK, zV are column vectors of monomials used to
%            specify K(x), V(x) in the vector form, V(x) = zV(x)'*C.
%            [default: zV = monomials(x,2); zK empty]
%
%       Monomials can alternatively specified for each subdomain of a
%       piecewise system:
%       -z1i:   Cell of monomials for shape multipliers
%       -z2i:   Cell of monomials for gradient multipliers
%       -zgi:   Cell of monomials for constraint multipliers
%       -zi:    Cell of monomials for boundary multipliers
%       -zKi:   Cell of monomials for multiple control functions
%       -zVi:   Cell of monomials for multiple Lyapunov functions
%
%       Advanced options:
%       -NstepBis:  Number of bisection steps [default = 1]
%       -betamax:   Maximum value of beta in beta-s1 step [default = 100]
%       -gammamax:  Maximum value of gamma in gamma-s2 step [default = 100]
%       -sosopt:    Options for feasibility steps [default = sosoptions]
%       -gsosopt:   Options for bisection steps [default = gsosoptions]
%                   See SOSOPTIONS and GSOSOPTIONS for more details.
%
%       Display, debug, and log
%       -display:   Display information of each iteration step 
%       -debug:     Display additional information during iteration
%                   [allowable values: 'on' or 'off'; default = 'off']
%       -log:       Save results to MAT files
%                   [allowable values: 'result' (store only final result),
%                   'step' (store results of each iteration step), or
%                   'none' (don't store results); default = 'none']
%       -logpath:   Specify path for result MAT files
%                   [default = 'data/']
%                   Note: the log path can be specified as character array
%                   or cell vector of character arrays; the first (or only)
%                   element is considered to be a folder. If this folder
%                   does not exist yet, it is created; if it does exist,
%                   any MAT files already stored inside is deleted.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2018-05-22
% * Changed:    2019-11-01
%
%% See also
%
% See ROAOPTIONS
%%

properties
    %f -- inherited from ROAOPTIONS
    %x -- inherited form ROAOPTIONS
    f1;
    f2;
    phi;
    xi;
    %z1 -- inherited from ROAOPTIONS
    %z2 -- inherited from ROAOPTIONS
    %zg -- inherited from CONROAOPTIONS
    zi;
    z1i;
    z2i;
    zgi;
    %zV  -- inherited from ROAOPTIONS
    %Vin -- inherited from ROAOPTIONS
    zVi;
    Vi0;
    %zK  -- inherited from CONROAOPTIONS
    %Kin -- inherited from CONROAOPTIONS
    zKi;
    Ki0;
    
    %display -- inherited from ROAOPTIONS
    debug = 'off';
    
    gammacheck = 'none';
end

methods
    function opt = pwroaoptions(f1, f2, varargin)
        % pwroaoptions(f, x, [u], ...)
        if ispvar(f2) && ...
                (nargin < 3 || ispvar(varargin{1}) || ischar(varargin{1}))
            % fall back to single ROA estimation
            x = f2;
            
            phi = -Inf;
            f2 = [];
        else
            phi = varargin{1};
            x = varargin{2};
            
            varargin(1:2) = [];
        end
            
        
        opt@conroaoptions(f1, x, varargin{:});
        
        opt.f1  = f1;
        opt.f2  = f2;
        opt.phi = phi;
        
        if isempty(opt.xi)
            opt.xi = opt.x;
        end
        
        if isempty(opt.zi)
            % XXX More intelligent selection?
            opt.zi =  monomials(opt.xi, 0:2);
        end
        
        if isempty(opt.zVi)
            opt.zVi = opt.zV;
        end
        
        if isempty(opt.z1i)
            opt.z1i = opt.z1;
        end
        
        if isempty(opt.z2i)
            opt.z2i = opt.z2;
        end
        
        if isempty(opt.zgi)
            opt.zgi = opt.zg;
        end
        
        if isempty(opt.zKi)
            opt.zKi = opt.zK;
        end
        
        if isempty(opt.Vi0) && ~isempty(opt.Vin)
            opt.Vi0 = {opt.Vin};
        end
        
        if isempty(opt.Ki0) && ~isempty(opt.Kin)
            opt.Ki0 = {opt.Kin};
        end
        
%         if isempty(opt.gammacheck) && length(opt.zVi) > 1
%             % default: activate gamma feasibility check
%             % for multiple Lyapunov functions
%             opt.gammacheck = 'feas';
%         elseif isempty(opt.gammacheck)
%             opt.gammacheck = 'none';
%         end

        if isempty(opt.logpath) && ~strcmp(opt.log,'none')
            opt.logpath = 'data/';
        end
    end
    
    % Set: zVi
    function opt = set.zVi(opt,value)
        if ismonom(value)
            opt.zVi = {value};
        elseif iscell(value) && ~isempty(value) && ismonom(value{1})
            opt.zVi = value;
        else
            error('zVi must be a non-empty cell of vectors of monomials.');
        end
        
        opt.zV = opt.zVi{1};
    end

    % Set: zKi
    function opt = set.zKi(opt,value)
        if isempty(value) || ismonom(value)
            opt.zKi = {value};
        elseif iscell(value) && ~isempty(value) && ismonom(value{1})
            opt.zKi = value;
        else
            error('zKi must be a non-empty cell of vectors of monomials.');
        end
        
        opt.zK = opt.zKi{1};
    end

    % Set: Vi0
    function opt = set.Vi0(opt,value)
        if iscell(value) && ~isempty(value) && isa(value{1},'polynomial')
            opt.Vi0 = value;
        else
            error('Lyapunov functions must be polynomials.');
        end
    end
    
    % Set: Ki0
    function opt = set.Ki0(opt,value)
        if iscell(value) && ~isempty(value) && isa(value{1},'polynomial')
            opt.Ki0 = value;
        else
            error('Controllers must be polynomials.');
        end
    end
    
    % Set: z1i
    function opt = set.z1i(opt,value)
        if ismonom(value) || isa(value,'double')
            opt.z1i = {value};
        elseif iscell(value) && ~isempty(value) && ismonom(value{1}) && ~isa(value{1},'double')
            opt.z1i = value;
        else
            error('Multiplier of piecewise beta step must be a monomial or constant.');
        end
        
        opt.z1 = opt.z1i{1};
    end
    
    % Set: z2i
    function opt = set.z2i(opt,value)
        if ismonom(value) && ~isa(value,'double')
            opt.z2i = {value};
        elseif iscell(value) && ~isempty(value) && ismonom(value{1}) && ~isa(value{1},'double')
            opt.z2i = value;
        else
            error('Multiplier of piecewise gamma step must be a monomial and non-constant.');
        end
        
        opt.z2 = opt.z2i{1};
    end 
    
    % Set: zgi
    function opt = set.zgi(opt,value)
        if isempty(value)
            opt.zgi = cell(1,1);
        elseif ismonom(value) && ~isa(value,'double')
            opt.zgi = {value};
        elseif iscell(value) && ~isempty(value) && ismonom(value{1}) && ~isa(value{1},'double')
            opt.zgi = value;
        else
            error('Multiplier of piecewise constraint step must be a monomial and non-constant.');
        end
        
        opt.zg = opt.zgi{1};
    end
    
    % Set: zi
    function opt = set.zi(opt,value)
        if ismonom(value) && ~isa(value,'double')
            opt.zi = {value};
        elseif iscell(value) && ~isempty(value) && ismonom(value{1}) && ~isa(value{1},'double')
            opt.zi = value;
        else
            error('Multiplier for boundary condition must be a monomial and non-constant.');
        end
    end
    
    % Set: xi
    function opt = set.xi(opt,value)
        if ispvar(value)
            opt.xi = value;
        else
            error('State vector for boundary condition must be polynomial variables.');
        end
    end
    
    % Set: gammacheck
    function opt = set.gammacheck(opt,value)
        AllowableVal = {'none' 'feas' 'check'};
        if ischar(value) && any(strcmp(value,AllowableVal))
            opt.gammacheck = value;
        else
            error('gammacheck must be one of ''none,'' ''feas,'' or ''check''.');
        end
    end
    
    % Set: debug
    function opt = set.debug(opt,value)
        switch value
            case 'off'
                % nothing to do
            case 'on'
                opt.display = 'on';
            otherwise
                error('debug can be ''on'' or ''off''. ');
        end
        
        % if not error
        opt.debug = value;
    end    
end

end