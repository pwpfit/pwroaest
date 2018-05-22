function [gbnds,sout1,sout2,info] = pwpcontain(pa,pb,p2,phi,z,zi,opts)
% Maximizes g subject to the piecewise set containment constraint
%   {x: p2(x) <= g} intersecting {phi(x) <= 0} is subset of {x: pa(x) <= 0}
%   {x: p2(x) <= g} intersecting {phi(x) >= 0} is subset of {x: pb(x) <= 0}
%
%% Usage & description
%
% Inputs:
%       -pa,pb: polynomials to describe piecewise superset
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
% * Changed:    2018-05-22
%
%% See also
%
% See PCONTAIN
%%

% polynomial variables
x = unique([pa.varname; pb.varname; p2.varname; phi.varname]);

% default options
if ~exist('opts','var') || isempty(opts)
    opts = gsosoptions;
end

%% Call GSOSOPT
% to solve:
%       max g such that
%       s, si in SOS
%       -( pa + (g-p2)*s' - phi*si' ) in SOS
%       -( pb + (g-p2)*s" + phi*si" ) in SOS
%
% GSOSOPT solves minimization of t := -g
%       min t such that
%       s, si in SOS
%       p2*s' - pa + t*s' + phi*si' in SOS
%       p2*s" - pb + t*s" - phi*si" in SOS
t   = pvar('g');
s1  = sosdecvar('c',z);
s2  = sosdecvar('c',z);
si1 = sosdecvar('c',zi);
si2 = sosdecvar('c',zi);

sosc(1) = s1 >= 0;
sosc(2) = s2 >= 0;
sosc(3) = si1 >= 0;
sosc(4) = si2 >= 0;

sosc(5) = (p2*s1 - pa) + t*s1 + phi*si1 >= 0;
sosc(6) = (p2*s2 - pb) + t*s2 - phi*si2 >= 0;

gopts = opts;
gopts.minobj = -opts.maxobj;
gopts.maxobj = -opts.minobj;
if strcmp(opts.display,'on')
    gopts.display = 'pcontain'; % Undocumented syntax for proper display
end

[info,dopt] = gsosopt(sosc,x,t,gopts);

if ~isempty(info.tbnds)
    gbnds = -info.tbnds([2 1]);
    
    sout1 = [
        subs(s1, dopt)
        subs(s2, dopt)
    ];
    sout2 = [
        subs(si1, dopt)
        subs(si2, dopt)
    ];
else
    gbnds = [];
    sout1 = [polynomial([]); polynomial([])];
    sout2 = [polynomial([]); polynomial([])];
end
