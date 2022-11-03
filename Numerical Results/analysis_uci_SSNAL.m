%**************************************************************************
% Computes clusterpath for n = 200 to test whether everything is working
% correctly. Compared to output of other algorithms in n200_test.R
%**************************************************************************

clear all;
addpath(strcat(pwd,'/Matlab/Solver.p/'));
addpath(strcat(pwd,'/Matlab/utils/'));
addpath(strcat(pwd,"/Matlab/utils_added"));

% Load results generated in R
[values, rownames, colnames] = load_uci_results("Output/uci_results.csv");

for i = 1:3
    pathX = strcat("Data/UCI/", rownames(i), "_X_umap.csv");
    pathG = strcat("Data/UCI/", rownames(i), "_lambdas.csv");
    
    % Get input for SSNAL
    X = csvread(pathX)';
    gammas = csvread(pathG)';
    k = values(i, find(colnames(1, :) == """k""", 1, "first"));
    phi = values(i, find(colnames(1, :) == """phi""", 1, "first"));
    
    % Run SSNAL
    [output_ssnal, cpu_time, W] = cvxc_ssnal(X, gammas, k, phi);
    
    % Get the elapsed time
    values(i, find(colnames(1, :) == """time_ssnal""", 1, "first")) = ...
        cpu_time;
    
    % Save results which are extended with SSNAL timings
    save_uci_results("Output/uci_results.csv", values, rownames, colnames);
end
