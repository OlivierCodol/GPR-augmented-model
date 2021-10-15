function h = plot_perf(h, xdata, ydata, zdata, zlabelstr, Q)

% sort input data
Z = reshape(zdata,[],length(xdata));

% filled contour plot
contourf( h , xdata , ydata , Z );
% imagesc( h , X , Y , Z ,'AlphaData',double(~isnan(Z)));
% h.YDir='normal';
ylim([min(ydata) max(ydata)])

% mark best Q value position
x = sort( repmat(xdata, 1, length(ydata)) );
y =       repmat(ydata, 1, length(xdata));
ix = find(Q==max(Q));
hold on
plot( x(ix) , y(ix) , 'k*','markersize',10,'Parent',h)

% add colorbar
colormap(h, 'parula')
c = colorbar(h);

% add variable name
pos_x = c.Position(1);
pos_y = c.Position(2) + c.Position(4) + 0.02*(c.Position(4)-c.Position(2));
annotation(h.Parent,'textbox',[pos_x ,pos_y,.05,.05],'string',zlabelstr,...
    'Interpreter','latex',...
    'String',zlabelstr,...
    'EdgeColor','none',...
    'VerticalAlignment','bottom',...
    'HorizontalAlignment','left',...
    'FontSize',10)




