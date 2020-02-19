function [beta,V,gamma,K,iter] = roaload(varargin)
% Loads stored ROA estimation results.
%
%% Usage & description
%
%   [beta,V,gamma,...] = roaload(filename,path)
%   [...] = roaload(filename,ropts)
%   [...] = roaload(iter,ropts)
%   [...] = roaload(ropts)
%
% Inputs:
%   - path:     file path (absolute/relative)
%   - filename: file in log folder (see ropts)
%   - iter:     iteration in log folder (see ropts)
%   - ropts:    options for ROA estimation; used to determine log folder.
%               See PWROAOPTIONS for defaults.
%
%   If no log file is specified, the final result is loaded.
%
% Outputs:
%   See PWROAEST for output variables.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:torbjoern.cunis@onera.fr>
% * Created:    2019-01-21
% * Changed:    2019-10-28
%
%% See also
%
% See PWROAOPTIONS, PWROAEST
%%

for i=1:nargin
    arg = varargin{i};
    
    if ~exist('filename','var') && ischar(arg) && ~exist('iter','var')
        filename = arg;
    elseif ~exist('path','var') && ischar(arg)
        path = arg;
    elseif ~exist('iter','var') && isnumeric(arg)
        iter = arg;
    elseif ~exist('ropts','var') && isa(arg,'roaoptions')
        ropts = arg;
    end
end

if ~exist('filename','var'),    filename = [];  end
if ~exist('path','var'),        path = [];      end
if ~exist('iter','var'),        iter = [];      end
if ~exist('ropts','var'),       ropts = [];     end


if ~isempty(path)
    if ~iscell(path)
        path = {path};
    end
elseif ~isempty(ropts)
    path = ropts.logpath;
elseif ~isempty(filename)
    path = '.';
else
    ropts = pwroaoptions(polynomial, polynomial, polynomial, pvar('x'));
    path = ropts.logpath;
end

if ~isempty(filename)
    [~,~,ext] = fileparts(filename);
    if isempty(ext)
        filename = [filename '.mat'];
    end
elseif ~isempty(iter)
    filename = sprintf('iter%d.mat',iter);
else
    filename = 'result.mat';
end

S = load([path{:} filename], '-mat');

if isfield(S, 'Ibar')
    beta  = S.Ibar;
else
    beta  = S.beta;
end
V     = S.V;
gamma = S.gamma;

if isfield(S, 'K')
    K = S.K;
else
    K = [];
end

if isfield(S, 'iter')
    iter = S.iter;
else
    iter = S;
end

end