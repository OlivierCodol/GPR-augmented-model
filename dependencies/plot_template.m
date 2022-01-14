function plot_template(ha, template)

x = repmat([0:0.1:1]', 11, 1);
y = sort(x);

% plot(x(template(:)==1), y(template(:)==1), '.', 'color', 'k'); hold on
% plot(x(template(:)==2), y(template(:)==2), 'o', 'color', 'k', 'MarkerSize', 4)
% plot(x(template(:)==3), y(template(:)==3), '.', 'color', 'k', 'MarkerSize', 15)

plot(x(template(:)==1), y(template(:)==1), '.', 'color', 'k'); hold on
plot(x(template(:)==2), y(template(:)==2), 'o', 'color', 'k', 'MarkerSize', 4)
plot(x(template(:)==5), y(template(:)==5), '.', 'color', 'k', 'MarkerSize', 15)
plot(x(template(:)==4), y(template(:)==4), '^', 'color', 'k', 'MarkerSize', 4, 'MarkerFaceColor', 'k')
plot(x(template(:)==3), y(template(:)==3), '+', 'color', 'k', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
plot(x(template(:)==6), y(template(:)==6), '*', 'color', 'k', 'MarkerSize', 4, 'MarkerFaceColor', 'k')

xlim([-0.1 1.1])
ylim([-0.1 1.1])
ha.XAxis.TickValues = 0:0.2:1;
ha.YAxis.TickValues = 0:0.2:1;
xlabel('salience (channel 1)')
ylabel('salience (channel 2)')
axis square