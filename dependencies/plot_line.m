function plot_line(h, r)

% plot line of best ratio
x_end = h.XLim(2);
y_end = h.YLim(2);

x2 = min( y_end*r , x_end );
y2 = min( y_end , x_end/r );

line( [0 x2] , [0 y2] ,'linestyle','--','color',[0 0 0],'parent',h,'linewidth',2)
