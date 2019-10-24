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
    z1i;
    z2i;
    %zV  -- inherited from ROAOPTIONS
    %Vin -- inherited from ROAOPTIONS
    zVi;
    Vi0;
    
    %display -- inherited from ROAOPTIONS
    debug = 'off';
    
    gammacheck = 'none';
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
        
        if isempty(opt.z1i)
            opt.z1i = opt.z1;
        end
        
        if isempty(opt.z2i)
            opt.z2i = opt.z2;
        end
        
        if isempty(opt.Vi0) && ~isempty(opt.Vin)
            opt.Vi0 = {opt.Vin};
        end
        
%         if isempty(opt.gammacheck) && length(opt.zVi) > 1
%             % default: activate gamma feasibility check
%             % for multiple Lyapunov functions
%             opt.gammacheck = 'feas';
%         elseif isempty(opt.gammacheck)
%             opt.gammacheck = 'none';
%         end
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
    
    % Set: z1i
    function opt = set.z1i(opt,value)
        if  ismonom(value) || isa(value,'double')
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
        if  ismonom(value) && ~isa(value,'double')
            opt.z2i = {value};
        elseif iscell(value) && ~isempty(value) && ismonom(value{1}) && ~isa(value{1},'double')
            opt.z2i = value;
        else
            error('Multiplier of piecewise gamma step must be a monomial and non-constant.');
        end
        
        opt.z2 = opt.z2i{1};
    end 
    
    % Set: zi
    function opt = set.zi(opt,value)
        if  ismonom(value) && ~isa(value,'double')
            opt.zi = {value};
        elseif iscell(value) && ~isempty(value) && ismonom(value{1}) && ~isa(value{1},'double')
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