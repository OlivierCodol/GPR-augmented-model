function c = redgreyblue(m,varargin)

if (nargin<1 || isempty(m))
    m = size(get(gcf,'colormap'),1);
end
if nargin==2
    CENTER = varargin{1};
else
    CENTER = .5;
end

m1 = floor(m*CENTER);
topval = 0.85; % brightness of the grey (if 1, then it becomes white)

r = reshape(linspace(0,topval,m1),[],1);
g = r;
r = [r; topval*ones(m-m1,1)];
g = [g; flipud(reshape(linspace(0,topval,m-m1),[],1))];
b = flipud(reshape(linspace(0,topval,m-m1),[],1));
b = [topval*ones(m1,1); b];

c = [r g b]; 

