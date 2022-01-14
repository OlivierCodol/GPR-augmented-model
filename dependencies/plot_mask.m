function h = plot_mask(h, xdata, ydata, zdata )
% plot outline of region for which zdata>0

%------------------------
%-- sort input data -----
%------------------------
Z = reshape(zdata,[],length(xdata));
mask = Z'>0;


%------------------------------
%-- define mask boundaries ----
%------------------------------
C = bwlabeln(mask);  % Recode each cluster with a specific number
Cv = C(:);          % reshape into a column vector
nC = unique(Cv(mask(:)==1)); % Count how many individual clusters

boundary_array = cell(numel(nC),1); % boundaries array

for k = 1:numel(nC)
    % need to add a line for the "diff" later.
    % if not, the row index would be deported by 1. Also, values on the 1st
    % row would not be detectable.
    Cn = [false(1,size(C,2)); C==nC(k)];
    [row,col]=find(diff(Cn)==1);
    if ~isempty(row)
        boundary_array{k} = bwtraceboundary(mask,[row(1), col(1)],'N');
    else
        boundary_array{k} = []; % I think that line is redundant...
    end
end

I = cellfun(@isempty,boundary_array); % Find empty cells
boundary_array = boundary_array(~I);               % Erase empty cells


%-----------------
%-- plot mask ----
%-----------------
for b = 1:numel(boundary_array)
    hold on
    mask_x = xdata( boundary_array{b}(:,1) );
    mask_y = ydata( boundary_array{b}(:,2) );
    plot( mask_x, mask_y, 'Color',[0 0 0],'LineWidth',1, 'Parent',h)
end
