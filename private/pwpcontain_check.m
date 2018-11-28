function [sout1,sout2,info] = pwpcontain_check(pa,p2,phi,z,zi,opts)
% Check piecewise set containment constraint
%   {x: p2(x) <= 0} intersecting (phi(x) <= 0} is subset of {x: pa(x) <= 0}
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
% * Created:    2018-11-28
% * Changed:    2018-11-28
%
%% See also
%
% See PWPCONTAIN
%%

% polynomial variables
x = unique([pa.varname; p2.varname; phi.varname]);

% default options
if ~exist('opts','var') || isempty(opts)
    opts = gsosoptions;
end


%% Call SOSOPT
% to solve:
%       s, si in SOS
%       -( pa - p2*s - phi*si ) in SOS
%
s1  = sosdecvar2('c',z);
si1 = sosdecvar2('ci',zi);

sosc(1) = s1 >= 0;
sosc(2) = si1 >= 0;

sosc(3) = (p2*s1 - pa) +  phi   *si1 >= 0;

sopts = opts;
if strcmp(opts.display,'on')
    sopts.display = 'pcontain'; % Undocumented syntax for proper display
end

[info,dopt] = sosopt(sosc,x,sopts);

if info.feas
    sout1 = subs(s1,  dopt);
    sout2 = subs(si1, dopt);
else
    sout1 = [];
    sout2 = [];
end
