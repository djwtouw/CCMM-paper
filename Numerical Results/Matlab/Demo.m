%%*************************************************************************
%% Downloaded from: https://blog.nus.edu.sg/mattohkc/softwares/convexclustering/
%% Demo: SSNAL algorithms for convex clustering
%% Copyright (c) 2021 by
%% Defeng Sun, Kim-Chuan Toh and Yancheng Yuan 
%% Citation:
%% [1] D.F. Sun, K.C. Toh, and Y.C. Yuan, Convex clustering: model, 
%%     theoretical guarantee and efficient algorithm, 
%%     Journal of Machine Learning Research, 22(9), 2021.
%% [2] Y.C. Yuan, D.F. Sun, and K.C. Toh, An efficient semismooth Newton 
%%     based algorithm for convex clustering, ICML 2018.
%%*************************************************************************
clear all;
data_dir = pwd;
addpath(strcat(pwd,'/Solver.p/'));
addpath(strcat(pwd,'/utils/'));
%% options for the demo
%%=========================================
%% Dataset NO.:
%% n is the number of data points
%% 1 Two Half Moon (default n = 1000)
%% 2 MNIST (n = 1000)
%% 3 unbalanced gaussian (n = 8500)
%% 4 Fisher_Iris (n = 150)
%% 5 Semisphere (n = 5000)
%% 6 Wine (n = 178)
%%==========================================
dataset_NO = 1;
if (dataset_NO==1)
   npts = 1000; 
   load(strcat(pwd,['/Data/Two_Half_Moon/Half_Moon_Balance_',int2str(npts),'.mat']));
elseif (dataset_NO==2)
   load(strcat(data_dir, '/Data/MNIST/mnist.mat'));
elseif (dataset_NO == 3)    
   load(strcat(data_dir, '/Data/Unbalanced_Gauss/unbalanced_gauss.mat'));
elseif (dataset_NO == 4)
   load(strcat(data_dir,'/Data/Fisher_Iris/Fisher_Iris.mat'));
elseif (dataset_NO == 5)
   load(strcat(data_dir,'/Data/Semisphere/semisphere.mat'));
elseif (dataset_NO == 6)
    load(strcat(data_dir, '/Data/Wine/wine_data.mat'));
end 
dataMatrix = X; %%d by (n=number of points)
%% normalize X
[dim.d,dim.n] = size(dataMatrix);
k_n = 10;
phi = 0.5;
%% Compute weights
fprintf('\n Start to compute weights.\n');
[weightVec1,NodeArcMatrix] = compute_weight(dataMatrix,k_n,phi,1);
%% Construct Amap
A0 = NodeArcMatrix;
Ainput.A = A0;
Ainput.Amap = @(x) x*A0;
Ainput.ATmap = @(x) x*A0';
Ainput.ATAmat = A0*A0'; %%graph Laplacian
Ainput.ATAmap = @(x) x*Ainput.ATAmat;
dim.E = length(weightVec1);
options.stoptol = 1e-6; %% tolerance for terminating the algorithm
options.num_k = k_n; %%number of nearest neighbors
%%==============================================
%% Stopping Criteria 
%% use_kkt == 1: based on relative KKT residual
%% use_kkt == 0: based on relative duality gap
%% Note that we will always terminate fast AMA based on relative duality gap
%%==============================================
options.use_kkt = 0;
%%===============================================
%% Implemented Algorithms
%% run_fastama == 1 : run fast AMA 
%% Reference: Eric C. Chi & Kenneth Lange, Splitting methods for convex clustering, JCGS 2015
%% run_admm == 1: run ADMM (used as Phase I of SSNAL if warmstart applied)
%% run_ssnal == 1: run SSNAL
run_fastama = 0;
run_admm  = 0;
run_ssnal = 1;

flag_save_solution = 0;
gamma_list = [10:-0.2:1];
if (run_fastama == 1)
    output_fastama = {};
    iter = 0;
    for gamma = gamma_list 
        iter = iter + 1;
        weightVec = gamma*weightVec1;
        %% max iteration number for fast AMA 
        options.maxiter = 20000;
        if (iter == 1)
            [Z_fastAMA, X_fastAMA, info_fastAMA] = ...
             fast_AMA(Ainput,dataMatrix,dim,weightVec,options);
        else
            [Z_fastAMA, X_fastAMA, info_fastAMA] = ...
             fast_AMA(Ainput,dataMatrix,dim,weightVec,options,Z_fastAMA);
        end
        if flag_save_solution == 1
            if exist('solution','var')
                clear 'solution';
            end
            solution.X = X_fastAMA;
            solution.Z = Z_fastAMA;
            info_fastAMA.solution = solution;
        end
        output_fastama{iter} = info_fastAMA;
    end
end

if (run_admm == 1)
    time_used_ADMM = 0;
    output_admm = {};
    iter = 0;
    options.use_kkt = 1;
    options.solver = 'direct';
    for gamma = gamma_list
        iter = iter + 1;
        weightVec = gamma*weightVec1;
        %% max iteration number for ADMM 
        options.maxiter = 10000;
        if (iter == 1)
            [obj_ADMM,Y_ADMM,X_ADMM,Z_ADMM,info_ADMM,runhist_ADMM] = ...
             ADMM(Ainput,dataMatrix,dim,weightVec,options);
            time_used_ADMM = time_used_ADMM + info_ADMM.time;
        else
            %% Warmstart by using the solution for the previous gamma 
            [obj_ADMM,Y_ADMM,X_ADMM,Z_ADMM,info_ADMM,runhist_ADMM] = ...
             ADMM(Ainput,dataMatrix,dim,weightVec,options,Y_ADMM,X_ADMM,Z_ADMM);
            time_used_ADMM = time_used_ADMM + info_ADMM.time;
        end
        if flag_save_solution == 1
            if exist('solution','var')
                clear 'solution';
            end
            solution.X = X_ADMM;
            solution.Y = Y_ADMM;
            solution.Z = Z_ADMM;
            info_ADMM.solution = solution;
        end
        output_admm{iter} = info_ADMM;
    end
end

if (run_ssnal == 1)
    time_used = 0;
    output_ssnal = {};
    iter = 0;
    options.use_kkt = 1;
    for gamma =  gamma_list
        iter = iter + 1;
        weightVec = gamma*weightVec1;
        %% max iteration number for SSNAL 
        options.maxiter = 100;
        if (iter == 1)
            %% Warmstart SSNAL with ADMM (optional)
            options.admm_iter = 50;
            [obj_SSNAL,Y_SSNAL,X_SSNAL,Z_SSNAL,info_SSNAL,runhist_SSNAL] = ...
             SSNAL(Ainput,dataMatrix,dim,weightVec,options);
            time_used = time_used + info_SSNAL.time;
        else
            %% Warmstart by using the solution for the previous gamma
            options.admm_iter = 20;
            [obj_SSNAL,Y_SSNAL,X_SSNAL,Z_SSNAL,info_SSNAL,runhist_SSNAL] = ...
             SSNAL(Ainput,dataMatrix,dim,weightVec,options,Y_SSNAL,X_SSNAL,Z_SSNAL);
            time_used = time_used + info_SSNAL.time;
        end
        tolClustering = 10*options.stoptol;
        [cluster_id, num_cluster] = find_cluster(X_SSNAL,tolClustering);
        info_SSNAL.NumCluster = num_cluster;
        info_SSNAL.ClusterId  = cluster_id; 
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

  