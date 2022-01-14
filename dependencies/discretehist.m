function h = discretehist(data, h)

[ count , valc ] = histcounts( categorical(data) );
val = cellfun(@str2num, valc);
x = [ val ; val ];
y = [ zeros(1, numel(count)) ; count ];

plot(x,y, 'color', 'k', 'parent', h); hold on
plot(val, count, '.', 'color', 'k','MarkerSize',12, 'parent', h)
plot([median(data) median(data)], [0 max(count)*2], 'LineStyle', '--', 'LineWidth', 2, 'parent', h, 'Color','r')

r = max(val) - min(val);
xlim([min(val) - r / 20  max(val) + r / 20])
ylim([0 max(count)*1.1])
h.XAxis.TickValues = val;
h.XAxis.TickLabels = cellfun(@(x) floor(str2double(x) * 1000) / 1000, valc, 'UniformOutput',false);
h.XAxis.TickLabelRotation = 25;

end
