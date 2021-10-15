# GPR-augmented-model
GRP model augmented with differential D1 and D2 functional weights.

This is the code and data for the article available here:
https://www.biorxiv.org/content/10.1101/2020.11.12.380451v1

## How to run the simulations to produce the data

Use ```run_simulations.m```.

That script runs the GPR model simulation for all combination of
parameter values using a grid approach. Consequently this is very
time-consuming. The simulations are already performed and the results are
saved in the ```data\``` folder. The simulation code is provided mainly to show how the
simulations were done, but feel free to run the simulations again if you want to.

Note: the full grid run took ~8h using a parfor loop.

## How to reproduce the analyses and figures from the article


Use ```data_analysis.m```.

That script runs all the analyses from the paper and plots all the
figures based on the data. Since the data is already included in this repository, you do NOT
need to run ```run_simulations.m``` before running ```data_analysis.m```.
