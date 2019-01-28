function [gbnds,sout1,sout2,info] = pwpcontain(pa,p2,phi,z,zi,opts)
% Maximizes g subject to the piecewise set containment constraint
%   {x: p2(x) <= g} intersecting {phi(x) <= 0} is subset of {x: pa(x) <= 0}
%
%% Usage & description
%
% Inputs:
%       -pa:    polynomials to describe piecewise superset
%       -p2:    polynomials to describe subset
%       -phi:   polynomial to describe boundary condition
%       -z,zi:  Nz-by-1 column vectors of monomials used to specify the 
%               SOS decision variables s(x), si(x) in the Gram matrix form, 
%               s(x)=z(x)'*C*z(x); si(x)=zi(x)'*C*zi(x).
%       -opts:  options for optimization; see GSOSOPTIONS.
%
% Outputs:
%       -gbnds: 1-by-2 vector [glb,gub] providing lower and upper bound on
%               the maximum value of g; will be empty if no feasible lower
%               bound is found.
%       -s,si:  2-by-1 vectors of multiplier functions proving piecewise
%               set containment with glb; will be empty if no feasible
%               lower bound is found.
%       -info:  information structure returned by GSOSOPT
%
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2018-05-22
% * Changed:    2018-05-23
%
%% See also
%
% See PCONTAIN
%%

% polynomial variables
x = unique([pa.varname; p2.varname; phi.varname]);

% default options
if ~exist('opts','var') || isempty(opts)
    opts = gsosoptions;
end


%% Call GSOSOPT
% to solve:
%       max g such that
%       s, si in SOS
%       -( pa + (g-p2)*s - phi*si ) in SOS
%
% GSOSOPT solves minimization of t := -g
%       min t such that
%       s, si in SOS
%       p2*s - pa + t*s + phi*si in SOS
t   = pvar('g');

if size(pa,1) == size(p2,1)
    s1  = sosdecvar2('c',z);
    si1 = sosdecvar2('ci',zi);
else
    s1  = sosmdecvar('c',z,size(pa,1)/size(p2,1));
    si1 = sosmdecvar('ci',zi,size(pa,1)/size(p2,1));
end

sosc = cell(3,1);

sosc{1} = s1 >= 0;
sosc{2} = si1 >= 0;

sosc{3} = (p2*s1 - pa) + t*s1 +  phi   *si1 >= 0;

sosc = vertcat(sosc{:});

gopts = opts;
gopts.minobj = -opts.maxobj;
gopts.maxobj = -opts.minobj;
if strcmp(opts.display,'on')
    gopts.display = 'pcontain'; % Undocumented syntax for proper display
end

[info,dopt] = gsosopt(sosc,x,t,gopts);

if ~isempty(info.tbnds)
    gbnds = -info.tbnds([2 1]);
    
    sout1 = subs(s1, dopt);
    sout2 = subs(si1, dopt);
else
    gbnds = [];
    sout1 = polynomial;
    sout2 = polynomial;
end
