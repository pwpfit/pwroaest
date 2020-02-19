classdef conroaoptions < droaoptions
% Options for ROA estimation under constraint control.
%
%% Usage & description
%
%   opts = conroaoptions(f, x, [u], ...)
%
% with inputs
%       -f: polynomial vector field
%       -x: state-space vector as PVAR
%       -u: input vector as PVAR
%
%   Name, Value:
%       -zK: Vector of monomials for control function;
%            i.e., zK is a column vector of monomials used to specifiy K(x)
%            in the vector form, K(x) = C'*zK(x).
%       -c:  constraint polynomial, vector field
%       -zg: Monomials for constraint multiplier; zg is a column vector of
%            monomials used to specifiy sg(x) in the Gram matrix form.
%            [default = monomials(x, 1:2)]
%       -Kin: Initial control function; default: [].
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2018-12-01
% * Changed:    2019-10-28
%
%% See also
%
% See ROAOPTIONS
%%

properties
    %f -- inherited from ROAOPTIONS
    %x -- inherited from ROAOPTIONS
    c;
    u;
    zg;
    
    zK;
    Kin;
end

methods
    function opt = conroaoptions(f, x, u, varargin)
        if nargin < 3
            u = [];
        elseif ischar(u)
            varargin = [{u} varargin];
            u = [];
        end
        
        opt@droaoptions(f, x, varargin{:});
        
        if isempty(opt.u)
            opt.u = u;
        end
        
        if ~isempty(u) && ~isempty(intersect(opt.x.varname, opt.u.varname))
            warning('x and u should be disjunct.');
        end
        %TODO: check f is linear in u if zK is given.

        % set zg to default only if constraint is given
        if isempty(opt.c)
            opt.c = polynomial(-1);
            %opt.zg = [];
        elseif isempty(opt.zg)
            opt.zg = monomials(opt.x, 0:2);
        end

    end
    
    function opt = set.c(opt,value)
        if isa(value,'polynomial')
            opt.c = value;
        else
            error('Constraint function must be polynomial.');
        end
    end
    
    function opt = set.zg(opt,value)
        if isempty(value)
            opt.zg = [];
        elseif ismonom(value)
            opt.zg = value;
        else
            error('zg must be vector of monomials.');
        end
    end
end

end