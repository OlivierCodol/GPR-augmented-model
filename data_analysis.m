% -----------------------------------------------
% This script runs all the analyses from the paper and plots all the
% figures based on the data. See 'run_simulations.m' for how the data was
% acquired.
% -----------------------------------------------
% O.Codol - codol.olivier@gmail.com
% 14-Oct-2021
% -----------------------------------------------



clearvars
close all
addpath('dependencies')
rootdir = getRootDirectory();
cd(rootdir)

% run_simulations


%% COMPUTE PERFORMANCE PROXIES

% simulated models values
load([rootdir '/data/augmented_models.mat'], 'hard', 'soft', 'Rw', 'Wd1', 'Wd2')
[Wx, d_Fh, d_Fs, H_max, S_max ] = compute_proxies(hard , soft, Rw);

% control model values
load([rootdir '/data/original_model.mat'],'hard_o','soft_o','Rw');
[Wx_o, d_Fh_o, d_Fs_o, H_max_o, S_max_o ] = compute_proxies(hard_o , soft_o, Rw);



%% COMPUTE RELATIVE PERFORMANCE & MERIT

% obtain ratio and tilde performance
[R_Wx  , T_Wx  ] = compute_perf( Wx , Wx_o );
[R_d_Fh, T_d_Fh] = compute_perf( d_Fh , d_Fh_o );
[R_d_Fs, T_d_Fs] = compute_perf( d_Fs , d_Fs_o );
[RH_max, TH_max] = compute_perf( H_max , H_max_o );
[RS_max, TS_max] = compute_perf( S_max , S_max_o );


% obtain merit Q
prod = R_Wx .* RH_max .* RS_max .* R_d_Fh .* R_d_Fs;
Q = log10(prod);
Q(prod==0) = nan;


%% ANALYZE RESULTS OF SIMULATIONS WITH STOCHASTIC WEIGHTS


load('data\noisy_models.mat', 'featureset')
F = featureset(~isnan(featureset(:,1)),:);

D1n = F(:, end-1);
D2n = F(:, end);
Rn = D1n ./ D2n;
Qn = F(:, end-2);

hf(1) = figure(1); clf
hf(1).Position(3) = 550;
hf(1).Position(4) = 200;


h = subplot(1,2,1);
discretehist(D1n, h);
xlabel('D1 weight')
ylabel('count')
legendline(hf(1), {'median'}, [1 0 0], 'position', [0.135 0.87 0.2 0.019], 'linestyle',{'--'}, 'linewidth',2);


h = subplot(1,2,2);
h = discretehist(D2n, h);
ix = [2:3:numel(h.XAxis.TickValues), 3:3:numel(h.XAxis.TickValues)];
h.XAxis.TickValues(ix) = [];
h.XAxis.TickLabels(ix) = [];
xlabel('D2 weight')
ylabel('count')
legendline(hf(1), {'median'}, [1 0 0], 'position', [0.58 0.87 0.2 0.019], 'linestyle',{'--'}, 'linewidth',2);



%% RATIO ANALYSIS

inv_W_D2 = 1 ./ Wd2';                       % inverse of w_d2
outer_product = inv_W_D2 * Wd1;             % all ratio combinations
ratio = reshape(outer_product, [],1);       % vectorise
Q_max = max(Q(:));                          % best Q value
ratio_max = ratio(Q(:)==Q_max);             % ratio for best Q value


% ratio plot
hf(2) = figure(2); clf;
hf(2).Units = 'pixels';
hf(2).Position = [800 20 600 300];
h = subplot(1,2,1);

semilogx( ratio, Q(:), '.','parent',h,'color','k')
xline( 1 ,'linestyle','--','color',[0 0 0],'parent',h)
yline( Q_max ,'linestyle','--','color',[0 0 0],'parent',h)
xlabel('W_{D1} / W_{D2}')
ylabel('Q value')
text(h.XLim(2)+20, Q_max,['Qmax'])
% text(15, Q_max+0.4,['Q = ' num2str(Q_max)])
% text(ratio_max+1, -1.1,      ['ratio = ' num2str(ratio_max)])
h.XTickLabel = {'0.1';'1';'10';'100';'1000'};
h.YAxis.Limits(2) = 0.5;

h = subplot(1,2,2);
histogram(D1n ./ D2n, 'edgecolor','none');
xline(1)
xline(median(D1n ./ D2n), '--', 'color','r', 'LineWidth',2)
ylim([0 max(histcounts(D1n ./ D2n)) * 1.1])
h.XAxis.TickDirection = 'out';
% h.YAxis.Scale = 'log';
% h.YAxis.TickValues = [1, 10, 100, 1000];
% h.YAxis.TickLabels = {'1';'10';'100';'1000'};
% h.XAxis.TickValues = 0:(1 + max(Rn(Rn<Inf)));
% h.XAxis.TickLabels = cellfun(@(x) num2str(x), num2cell(h.XAxis.TickValues), 'UniformOutput',false);
xlabel('D1 / D2 weight ratio')
ylabel('count')
legendline(hf(2), {'median'}, [1 0 0], 'position', [0.75 0.87 0.2 0.019], 'linestyle',{'--'}, 'linewidth',2);


annotation(hf(2),'textbox',[.05 .97 .05 .05],'String','a','EdgeColor','none','FontWeight','bold','FontSize',14)
annotation(hf(2),'textbox',[.48 .97 .05 .05],'String','b','EdgeColor','none','FontWeight','bold','FontSize',14)



%% PLOT SURFACE PLOTS

hf(3) = figure(3); clf;
hf(3).Units = 'pixels';
hf(3).Position = [20 20 833 539];

h(1)=subplot(3,2,1); plot_perf(h(1), Wd1, Wd2, T_d_Fh, '$\Delta \tilde{F}_h$', Q);
h(2)=subplot(3,2,2); plot_perf(h(2), Wd1, Wd2, T_d_Fs, '$\Delta \tilde{F}_s$', Q);
h(3)=subplot(3,2,3); plot_perf(h(3), Wd1, Wd2, TH_max, '$\tilde{H}_{max}$'   , Q);
h(4)=subplot(3,2,4); plot_perf(h(4), Wd1, Wd2, TS_max, '$\tilde{S}_{max}$'   , Q);
h(5)=subplot(3,2,5); plot_perf(h(5), Wd1, Wd2, T_Wx  , '$\tilde{w}_x$'       , Q);
h(6)=subplot(3,2,6); plot_perf(h(6), Wd1, Wd2, Q     , 'Q'                   , Q);
                     plot_mask(h(6), Wd1, Wd2, Q     ); % outline region with Q>0

ylabel('W_{D2}','Rotation',0,'Parent',h(1))
ylabel('W_{D2}','Rotation',0,'Parent',h(3))
ylabel('W_{D2}','Rotation',0,'Parent',h(5))
xlabel('W_{D1}','Parent',h(5))
xlabel('W_{D1}','Parent',h(6))

annotation(hf(3),'textbox',[.05 .9 .05 .05],'String','a','EdgeColor','none','FontWeight','bold','FontSize',14)
annotation(hf(3),'textbox',[.05 .6 .05 .05],'String','c','EdgeColor','none','FontWeight','bold','FontSize',14)
annotation(hf(3),'textbox',[.05 .3 .05 .05],'String','e','EdgeColor','none','FontWeight','bold','FontSize',14)
annotation(hf(3),'textbox',[.51 .9 .05 .05],'String','b','EdgeColor','none','FontWeight','bold','FontSize',14)
annotation(hf(3),'textbox',[.51 .6 .05 .05],'String','d','EdgeColor','none','FontWeight','bold','FontSize',14)
annotation(hf(3),'textbox',[.51 .3 .05 .05],'String','f','EdgeColor','none','FontWeight','bold','FontSize',14)



%% PLOT EXAMPLE TEMPLATE FIT FUNCTIONS


% vectorise w_d1 and w_d2
Wd1_vec = sort( repmat(Wd1, 1, length(Wd2)) );
Wd2_vec =       repmat(Wd2, 1, length(Wd1));

% printing templates
fprintf(['\n\tmodel:\t\tH_max:'...
                   '\t\tS_max:'...
                   '\t\td_Fh:'...
                   '\t\td_Fs:'...
                   '\t\tWx:'...
                   '\t\t\tQ:'...
                   '\t\t\tWd1:'...
                   '\t\tWd2:\n'])
formatSpec    = '\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n';
formatSpecNaN = '\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t\t%.2f\t\t\t%.2f\t\t\t%.2f\t\t%.2f\n';


               
% create figure
hf(4) = figure(4); clf;
hf(4).Units = 'pixels';
hf(4).Position = [800 420 412 385];


%--------------------------------
%-- control model ----
subplot(2,2,1)
plot(Rw, hard_o,'color','k','LineStyle','-' ); hold on
plot(Rw, soft_o,'color','k','LineStyle','-.');
text(7, 8,'Q = 0.00','FontSize',8)
grid minor
ylabel('fit to template (%)')
xlim([1 10])
ylim([0 100])

fprintf('\torig')
fprintf(formatSpec, H_max_o , S_max_o , d_Fh_o , d_Fs_o , Wx_o , 0 , 1 , 1 )


%--------------------------------
%-- failed model ----
subplot(2,2,2)
ix=1151;
plot(Rw, hard(:,ix),'color','k','LineStyle','-' ); hold on
plot(Rw, soft(:,ix),'color','k','LineStyle','-.');
text(5.5, 8,'Q = undefined','FontSize',8)
grid minor
xlim([1 10])
ylim([0 100])

fprintf('\tfail')
fprintf(formatSpecNaN, H_max(ix) , S_max(ix) , d_Fh(ix)    , d_Fs(ix) ,...
                       Wx(ix)    , Q(ix)     , Wd1_vec(ix) , Wd2_vec(ix) )


%--------------------------------
%-- best model ----
subplot(2,2,3)
[~,ix] = max(Q);
plot(Rw, hard(:,ix),'color','k','LineStyle','-' ); hold on
plot(Rw, soft(:,ix),'color','k','LineStyle','-.');
text(5.8, 8,['Q = ' num2str(Q(ix))],'FontSize',8)
grid minor
ylabel('fit to template (%)')
xlabel('R_w')
xlim([1 10])
ylim([0 100])

fprintf('\tbest')
fprintf(formatSpec, H_max(ix) , S_max(ix) , d_Fh(ix)    , d_Fs(ix) ,...
                    Wx(ix)    , Q(ix)     , Wd1_vec(ix) , Wd2_vec(ix) )
                
                
%--------------------------------
%-- worse-than-control model ----
subplot(2,2,4)
ix=1104;
plot(Rw, hard(:,ix),'color','k','LineStyle','-' ); hold on
plot(Rw, soft(:,ix),'color','k','LineStyle','-.');
text(5.8, 8,['Q = ' num2str(Q(ix))],'FontSize',8)
grid minor
xlabel('R_w')
xlim([1 10])
ylim([0 100])

fprintf('\tworse')
fprintf(formatSpec, H_max(ix) , S_max(ix) , d_Fh(ix)    , d_Fs(ix) ,...
                    Wx(ix)    , Q(ix)     , Wd1_vec(ix) , Wd2_vec(ix) )

                
                
                
                
% plot legend & annotate panels
legendline(hf(4),{'soft mode';'hard mode'},...  % labels
                 zeros(2,3),...                 % colors
                 'linestyle',[{'-.'};{'-'}],... % linestyle
                 'position',[.26 .76 .3 .035]... % position in figure
             );
annotation(hf(4),'textbox',[.00 .95 .05 .05],'String','a','EdgeColor','none','FontWeight','bold','FontSize',12)
annotation(hf(4),'textbox',[.47 .95 .05 .05],'String','b','EdgeColor','none','FontWeight','bold','FontSize',12)
annotation(hf(4),'textbox',[.00 .45 .05 .05],'String','c','EdgeColor','none','FontWeight','bold','FontSize',12)
annotation(hf(4),'textbox',[.47 .45 .05 .05],'String','d','EdgeColor','none','FontWeight','bold','FontSize',12)





%% PLOT TEMPLATES

%--------------------------------
%-- control model ----
hf(5) = figure(5); clf
hf(5).Position = [50 50 500 500];

[OPTI_HARD, OPTI_SOFT] = get_optimal_templates();
hard_soft_array = [ lambda2Rw(0.294) , lambda2Rw(0.818) ];
[~, ~, ~, out] = gpr_model(1, 1, hard_soft_array);


plot_template( subplot(2,2,1), out{1} )
plot_template( subplot(2,2,2), out{2} )
plot_template( subplot(2,2,3), OPTI_HARD )
plot_template( subplot(2,2,4), OPTI_SOFT )


annotation(hf(5),'textbox',[.05 .95 .05 .05],'String','a','EdgeColor','none','FontWeight','bold','FontSize',12)
annotation(hf(5),'textbox',[.47 .95 .05 .05],'String','b','EdgeColor','none','FontWeight','bold','FontSize',12)
annotation(hf(5),'textbox',[.05 .45 .05 .05],'String','c','EdgeColor','none','FontWeight','bold','FontSize',12)
annotation(hf(5),'textbox',[.47 .45 .05 .05],'String','d','EdgeColor','none','FontWeight','bold','FontSize',12)



%% PLOT GPi OUTPUT DYNAMICS

hf(6) = figure(6); clf
hf(6).Position = [50 50 552 250];

[salience, GPi_output] = dynamics_sim(1, 1, 2);
t = 0.01:0.01:3;

subplot(1,2,1)
plot(t, GPi_output(1,:), 'color', 'k'); hold on
plot(t, salience(1,:), '--', 'color', 'k')
ylim([-0.05 0.8])
ylabel('signal level (a.u.)')
xlabel('timestep (a.u.)')
title('Channel 1')

subplot(1,2,2)
plot(t, GPi_output(2,:), 'color', 'k'); hold on
plot(t, salience(2,:), '--', 'color', 'k')
ylim([-0.05 0.8])
xlabel('timestep (a.u.)')
title('Channel 2')


% plot legend
legendline(hf(6),{'GPi output';'';'';'salience'},...  % labels
                 zeros(4,3),...                 % colors
                 'linestyle',[{'-'};{'none'};{'none'};{'--'}],... % linestyle
                 'position',[.14 .8 .3 .075]... % position in figure
             );





%% Format, save and close all


thedir = [pwd '\figures\gpr_analysis_' date_n_time '\'];
mkthedir(thedir)
addpath('dependencies\altmany-export_fig-9676767')

for k = 1:numel(hf)
    f = hf(k);
    export_fig(f,[thedir num2str(k)], '-png', '-transparent', '-m2');
%     fig2svg([thedir num2str(k) '.svg'],f);
end

rmpath('dependencies\altmany-export_fig-9676767')
rmpath('dependencies')





