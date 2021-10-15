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
[hard_original, soft_original] = gpr_model(1, 1, Rw);
save([datadir '/original_model.mat'], 'hard_original', 'soft_original', 'Rw')



%-------------------------------
% augmented model
%-------------------------------
[hard, soft] = gpr_model(Wd1, Wd2, Rw);
save([datadir '/augmented_models.mat'], 'hard', 'soft', 'Rw', 'Wd1', 'Wd2')


