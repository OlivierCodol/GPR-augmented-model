function mkthedir(thedir,varargin)
if nargin>1
    warnme = varargin{1};
else
    warnme = 0;
end
if exist(thedir,'dir')~=7
    mkdir(thedir)
elseif warnme==1
    warning('The directory already exists!')
end
end