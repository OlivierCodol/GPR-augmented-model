% -----------------------------------------------
% This script runs the GPR model simulation for all combination of
% parameter values using a grid approach. Consequently this is very
% time-consuming. The simulations are already performed and the results are
% saved in the 'data' folder. This code is provided mainly to show how the
% simulations were done, but feel free to run the simulations again if you
% want to.
% See 'data_analysis.m' for how the data was analysed and the figures
% plotted.
% Note: The full grid run took ~8h using a parfor loop.
% -----------------------------------------------
% O.Codol - codol.olivier@gmail.com
% 14-Oct-2021
% -----------------------------------------------


clearvars
addpath('dependencies')
rootdir = getRootDirectory();
cd(rootdir)
datadir = [rootdir '/data'];
mkthedir(datadir)


%-------------------------------
% variables of interest
%-------------------------------
Wd1 = linspace(0,10,401);   % values D1 weigths will take
Wd2 = linspace(0,11/9,50);  % values D2 weigths will take
Rw = linspace(1,10,1000);   % weight ratio for tonic dopamine levels



%-------------------------------
% control model
%-------------------------------
[hard_o, soft_o] = gpr_model(1, 1, Rw);
save([datadir '/original_model.mat'], 'hard_o', 'soft_o', 'Rw')



%-------------------------------
% augmented model
%-------------------------------
[hard, soft] = gpr_model(Wd1, Wd2, Rw);
save([datadir '/augmented_models.mat'], 'hard', 'soft', 'Rw', 'Wd1', 'Wd2')
 


%-------------------------------
% noisy models
%-------------------------------
[Wx_o, d_Fh_o, d_Fs_o, H_max_o, S_max_o] = compute_proxies(hard_o , soft_o, Rw);

Rw_n = linspace(1,10,100);   % weight ratio for tonic dopamine levels
Wd1_n = linspace(0,3,20);   % values D1 weigths will take
Wd2_n = linspace(0,11/9,20);  % values D2 weigths will take

n_sim = 1000;
featureset = nan(n_sim, 8);
[Wd1_vec, Wd2_vec] = meshgrid(Wd1_n, Wd2_n);
W = struct;

for k = 1:n_sim
    fprintf(['  Model #' num2str(k) ' : '])
    [hard, soft, W_k] = gpr_model(Wd1_n, Wd2_n, Rw_n, 'weight_bounds', [0.95 1.05]);
    
    [Wx, d_Fh, d_Fs, H_max, S_max] = compute_proxies(hard , soft, Rw_n);
    
    R_Wx   = compute_perf( Wx , Wx_o );
    R_d_Fh = compute_perf( d_Fh , d_Fh_o );
    R_d_Fs = compute_perf( d_Fs , d_Fs_o );
    RH_max = compute_perf( H_max , H_max_o );
    RS_max = compute_perf( S_max , S_max_o );
    
    prod = R_Wx .* RH_max .* RS_max .* R_d_Fh .* R_d_Fs;
    Q = log10(prod);
    Q(prod==0) = nan;
    
    [~,ix] = max(Q);
    featureset(k,:) = [H_max(ix), S_max(ix), d_Fh(ix), d_Fs(ix), Wx(ix), Q(ix), Wd1_vec(ix), Wd2_vec(ix)];
    
    D1 = featureset(1:k, end-1);
    D2 = featureset(1:k, end);
    figure(2); clf
    subplot(1,3,1); histogram(D1)
    subplot(1,3,2); histogram(D2)
    subplot(1,3,3); histogram(D1./D2)
    drawnow;
    
    W = addstr(W,W_k, 'column');
    save([datadir '/noisy_models.mat'], 'featureset', 'W')
end



