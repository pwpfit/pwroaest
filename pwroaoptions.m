classdef pwroaoptions < roaoptions
% Options for piecewise ROA estimation.
%
%% Usage & description
%
%   opts = pwroaoptions(f1, f2, phi, x, ...)
%
% with inputs
%       -f1:  first polynomial vector field
%       -f2:  second polynomial vector field
%       -x:   state-space vector as PVAR
%       -phi: boundary condition (scalar field)
%
%   Name, Value:
%       -xi: variables of boundary condition as PVAR; must be subset of the
%            state-space variables. [default =  opts.x]
%       -zi: Monomials for boundary multiplier in gamma-s2-si step of V-s 
%            iteration. Specifically, zi is a column vector of monomials 
%            used to specify si(xi) in the Gram matrix form, 
%            si(xi)=zi(xi)'*C*zi(xi). [default =  monomials(xi, 1:2) ]
%       -zVi:   Cell of monomials for multiple Lyapunov functions.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2018-05-22
% * Changed:    2018-09-09
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
    zi;
    %zV  -- inherited from ROAOPTIONS
    %Vin -- inherited from ROAOPTIONS
    zVi;
    Vi0;
end

methods
    function opt = pwroaoptions(f1, f2, phi, x, varargin)
        opt@roaoptions(f1, x, varargin{:});
        
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
        
        if isempty(opt.Vi0) && ~isempty(opt.Vin)
            opt.Vi0 = {opt.Vin};
        end
    end
    
    % Set: zV
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


    % Set: Vin
    function opt = set.Vi0(opt,value)
        if iscell(value) && ~isempty(value) && isa(value{1},'polynomial')
            opt.Vi0 = value;
        else
            error('Lyapunov functions must be polynomials.');
        end
    end
    
    % Set: zi
    function opt = set.zi(opt,value)
        if  ismonom(value) && ~isa(value,'double')
            opt.zi = value;
        else
            error('Multiplier for boundary condition must be a monomial and non-constant.');
        end
    end
    
    % Set: xi
    function opt = set.xi(opt,value)
        if  ispvar(value)
            opt.xi = value;
        else
            error('State vector for boundary condition must be polynomial variables.');
        end
    end
end

end