%**************************************************************************
% Measures elapsed times for the computation of the clusterpath using the
% SSNAL algorithm.
%**************************************************************************

clear all;
addpath(strcat(pwd,'/Matlab/Solver.p/'));
addpath(strcat(pwd,'/Matlab/utils/'));
addpath(strcat(pwd,"/Matlab/utils_added"));


% Settings for hyperparameters
gammas = (0.0:0.2:110.0);
phi = 2.0;
k = 15;

% Select the sizes for the data
N = (1000:1000:5000);

% Create vector to store computation times
timings = zeros(length(N), 1);

% The base path to the folder where the data is stored
data_path_base = "Data/Half Moons/Half_Moons_";


% Iterate over the different sizes of the data sets
for n_i = (1:length(N))
    data_path_n_i = strcat(data_path_base, int2str(N(n_i)), "/X_");
    
    % Iterate over the different realizations of the data for each size
    for d_i = 1:10
        disp(strcat("Analyzing n = ", int2str(N(n_i)), ", d = ", int2str(d_i)))
        % Read in the data
        data_path = strcat(data_path_n_i, int2str(d_i), ".csv");
        X = csvread(data_path)';
        
        % Compute results SSNAL
        [output_ssnal, cpu_time, W] = cvxc_ssnal(X, gammas, k, phi);
        timings(n_i) = timings(n_i) + cpu_time;
        
        if N(n_i) == 1000
            clusterpath = transform_output(output_ssnal);
            
            % Write output
            output_path_base = "Output/SSNAL Clusterpaths/n1000_";
            output_path = strcat(output_path_base, int2str(d_i), ...
                "_ssnal_clusterpath.csv");
            writematrix(clusterpath, output_path);
        end
        
        if N(n_i) == 5000
            clusterpath = transform_output(output_ssnal);
            
            % Write output
            output_path_base = "Output/SSNAL Clusterpaths/n5000_";
            output_path = strcat(output_path_base, int2str(d_i), ...
                "_ssnal_clusterpath.csv");
            writematrix(clusterpath, output_path);
        end
        
        % Write the times to csv
        writematrix(timings, "Output/n1000_to_5000_timings_ssnal.csv")
    end
end