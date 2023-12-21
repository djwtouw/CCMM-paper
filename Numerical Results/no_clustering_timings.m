%**************************************************************************
% Measures elapsed times for the computation of the clusterpath using the
% SSNAL algorithm. Data of interest are sets ranging from 500 to 20000.
% Additionally, the attained objective value is used later by CCMM to
% determine whether it has converged.
%**************************************************************************

clear all;
addpath(strcat(pwd,'/Matlab/Solver.p/'));
addpath(strcat(pwd,'/Matlab/utils/'));
addpath(strcat(pwd,"/Matlab/utils_added"));


% Settings for hyperparameters
phi = 0.5;
k = 15;

% Select the sizes for the data
N = [500:500:20000];

% Number of repetitions
n_reps = 10;

% Number of lambdas
n_lambdas = 5;

% Create vector to store computation times
timings = zeros(length(N), 1);

% Iterate over the different sizes of the data sets
for n_i = (1:length(N))
    for r = 1:n_reps
        disp(strcat("Analyzing n = ", int2str(N(n_i)), ", r = ", int2str(r)));

        % Vector holding losses
        losses = zeros(n_lambdas, 1);

        % Read in the data
        fname = strcat("Data/Grid/X_", int2str(N(n_i)), "_", int2str(r), ".csv");
        X = csvread(fname)';

        % Read in lambda
        fname = strcat("Data/Grid/lambda_", int2str(N(n_i)), "_", int2str(r), ".bin");
        gammas = load_bin(fname)';

        % Compute results SSNAL
        [output_ssnal, cpu_time, W] = cvxc_ssnal(X, gammas, k, phi);
        timings(n_i) = timings(n_i) + cpu_time;

        % Collect objective values
        for i = 1:length(gammas)
            losses(i) = output_ssnal{1, i}.obj(1);
        end

        % Write the losses to a binary file
        fname = strcat("Output/Grid/losses_", int2str(N(n_i)), "_", int2str(r), ".bin");
        write_bin(losses, fname);
    end
    
    % Print timings
    display(timings);

    % Write the times to csv
    writematrix(timings, "Output/no_clustering_timings_ssnal.csv");
end
