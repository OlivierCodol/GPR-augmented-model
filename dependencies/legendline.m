function varargout = legendline(HF,LABELS,CC,varargin)
% Create a line legend, returns an axis handle
% Syntax is 
%           legendline(HF,LABELS,CC,varargin)
% HF is the figure handle
% LABELS is a cell array wih N labels
% CC is a 1*3 or N*3 rgb color matrix
%
% Options are:
%   'linewidth'     Scalar or N-element vector indicating width of lines
%                   (default 1).
%   'linestyle'     Style of each line, 1- or N-elements cell array.
%   'fontsize'      Size of text, scalar or N-elements vector.
%   'position'      4-elements position in figure (default [0 0.05 1 .95])
%   'outerposition' Outer position in figure
%   'fliplr'        set the lines on the left of the text (0, default) or
%                   on the right of the text (1).
%
%
% O.Codol 02 Feb. 2020
% codol.olivier@gmail.com
%---------------------------------------------------------



if ~iscell(LABELS); error('labels should be a 1*n or n*1 cell array'); end
LABELS = LABELS(:);

L       = parsevarargin(varargin,'linewidth',1);
LS      = parsevarargin(varargin,'linestyle',{'-'});
N       = size(LABELS,1);
P       = parsevarargin(varargin,'Position',[0 0.05 1 .95]);
FS      = parsevarargin(varargin,'fontsize',floor(500*P(4)/N));
OP      = parsevarargin(varargin,'OuterPosition',[]);
FLR     = parsevarargin(varargin,'fliplr',0);

if size(CC,1)<N;    CC = repmat(CC,   N,1);     end
if numel(L) <N;     L  = repmat(L(:), N,1);     end
if numel(LS)<N;     LS = repmat(LS(:),N,1);     end
if numel(FS)<N;     FS = repmat(FS(:),N,1);     end

h = axes('Parent',HF,'Position',P); axes(h);
if ~isempty(OP); h.OuterPosition = OP; end
h.Color='none';
h.XAxis.Visible='off';
h.YAxis.Visible='off';


% Flip sides or not
if FLR==1               % labels on the left
    xt=.81;
    X=[.85 .99];
    POS = 'right';
elseif FLR==0           % lines on the left
    xt=.19;
    X=[.01 .15];
    POS = 'left';
end


% For each label
for k=1:N
   C = CC(k,:);
   Y = ([.5 .5]+k-1)/N;
   line(X,Y,'parent',h,'color',C,'linewidth',L(k),'linestyle',LS{k});
   
   yt=(k-.45)/N;
   text(xt,yt,LABELS{k},'FontSize',FS(k),'HorizontalAlignment',POS)
end


xlim([0 1])

varargout = {h};

