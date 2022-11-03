%**************************************************************************
% Computes clusterpath for n = 200 to test whether everything is working
% correctly. Compared to output of other algorithms in 
% [Numerical Results/n200_test_CCMM_AMA.R]
%**************************************************************************

clear all;
addpath(strcat(pwd,'/Matlab/Solver.p/'));
addpath(strcat(pwd,'/Matlab/utils/'));
addpath(strcat(pwd,"/Matlab/utils_added"));

% Settings for hyperparameters
gammas = (0.0:0.2:10.0);
phi = 2.0;
k = 10;

% Load data
path = strcat("Data/Half Moons/Half_Moons_0200/X.csv");
X = csvread(path)';

% Run SSNAL
[output_ssnal, cpu_time, W] = cvxc_ssnal(X, gammas, k, phi);

% Transform the SSNAL output into a clusterpath of the CCMM form
clusterpath = transform_output(output_ssnal);

%% Write output
output_path = "Output/SSNAL Clusterpaths/n200_ssnal_clusterpath.csv";
writematrix(clusterpath, output_path);

