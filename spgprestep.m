function [gbnds,s] = spgprestep(V,f,H,x,z,zi,L2,opts)

% default options
if ~exist('opts','var') || isempty(opts)
    opts = gsosoptions;
end

% number of cells
N = length(f);

sosc = cell(1,N);

%% Call GSOSOPT
% to solve:
%       max g such that
%       s, si in SOS
%       -( pa + (g-p2)*s - H'*si ) in SOS
%
% GSOSOPT solves minimization of t := -g
%       min t such that
%       s, si in SOS
%       p2*s - pa + t*s + H'*si in SOS
t   = pvar('g');

s = cell(N,2);

for j=1:N
    % ensure H is column vector
    assert(iscolumn(H{j}), 'H must be row or column of polynomials');

    % number of constraints Hi
    if isempty(H{j})
        % fall back to polynomial containment problem
        % See PCONTAIN.
        k  = 0;
    else
        k = length(H{j});
    end

    % origin in cell
    if iscell(zi) && all(double(subs(H{j}, x, zeros(size(x)))) <= 0)
        zj = zi{1};
    elseif iscell(zi)
        zj = zi{end};
    end

    s0 = sosdecvar(['c' num2str(j)],z);
    sj = polynomial(zeros(size(H{j})));
    for i=1:k
        sj(i) = sosdecvar(['c' num2str(j) '_' num2str(i)],zj);
    end

    sosc{j} = [
        s0 >= 0
        sj >= 0
        (V*s0 - (jacobian(V,x)*f{j}+L2)) + t*s0 + nonempty(H{j}'*sj) >= 0
    ];

    s{j,1} = s0;
    s{j,2} = sj;
end

gopts = opts;
gopts.minobj = -opts.maxobj;
gopts.maxobj = -opts.minobj;
if strcmp(opts.display,'on')
    gopts.display = 'pcontain'; % Undocumented syntax for proper display
end

[info,dopt] = gsosopt(vertcat(sosc{:}),x,t,gopts);

if ~isempty(info.tbnds)
    gbnds = -info.tbnds([2 1]);
    
    for j=1:N
        s{j,1} = subs(s{j,1}, dopt);
        s{j,2} = subs(s{j,2}, dopt);
    end
else
    gbnds = [];
    s(:)  = {polynomial};
end

end

function a = nonempty(a)
    if isempty(a)
        a = 0;
    end
end