%**************************************************************************
% Function to perform convex clustering using the SSNAL algorithm. Based on
% the contents of Demo.m
%**************************************************************************


function [output_ssnal, time_used, weightVec1] = cvxc_ssnal(X, gamma_list, k_n, phi)
    [dim.d,dim.n] = size(X);

    % Compute weights
    fprintf('\n Start to compute weights.\n');
    [weightVec1, NodeArcMatrix] = compute_weight(X, k_n, phi, 1);

    % Construct Amap
    A0 = NodeArcMatrix;
    Ainput.A = A0;
    Ainput.Amap = @(x) x*A0;
    Ainput.ATmap = @(x) x*A0';
    Ainput.ATAmat = A0*A0';     %%graph Laplacian
    Ainput.ATAmap = @(x) x*Ainput.ATAmat;
    dim.E = length(weightVec1);
    options.stoptol = 1e-6;     %% tolerance for terminating the algorithm
    options.num_k = k_n;        %%number of nearest neighbors
    
    % Prevent printing for SSNAL
    options.printyes = 0;
    options.printminoryes = 0;

    flag_save_solution = 1;

    time_used = 0;
    output_ssnal = {};
    iter = 0;
    options.use_kkt = 1;
    for gamma =  gamma_list
        iter = iter + 1;
        weightVec = gamma*weightVec1;
        % max iteration number for SSNAL 
        options.maxiter = 100;
        if (iter == 1)
            % Warmstart SSNAL with ADMM (optional)
            options.admm_iter = 50;
            [obj_SSNAL,Y_SSNAL,X_SSNAL,Z_SSNAL,info_SSNAL,runhist_SSNAL] = ...
             SSNAL(Ainput,X,dim,weightVec,options);
            time_used = time_used + info_SSNAL.time;
        else
            % Warmstart by using the solution for the previous gamma
            options.admm_iter = 20;
            [obj_SSNAL,Y_SSNAL,X_SSNAL,Z_SSNAL,info_SSNAL,runhist_SSNAL] = ...
             SSNAL(Ainput,X,dim,weightVec,options,Y_SSNAL,X_SSNAL,Z_SSNAL);
            time_used = time_used + info_SSNAL.time;
        end
        tolClustering = 10*options.stoptol;
        [cluster_id, num_cluster] = find_cluster(X_SSNAL,tolClustering);
        info_SSNAL.NumCluster = num_cluster;
        info_SSNAL.ClusterId  = cluster_id;
        info_SSNAL.gamma = gamma;
        if (flag_save_solution == 1)
            if exist('solution','var'); clear 'solution'; end
            solution.X = X_SSNAL;
            solution.Y = Y_SSNAL;
            solution.Z = Z_SSNAL;
            info_SSNAL.solution = solution;
        end
        output_ssnal{iter} = info_SSNAL;      
    end
end
