%%********************************************************************
%% Compute weights for convex clustering.
%% Copyright (c) 2021 by
%% Defeng Sun, Kim-Chuan Toh and Yancheng Yuan 
%% Citation:
%% [1] D.F. Sun, K.C. Toh, and Y.C. Yuan, Convex clustering: model, theoretical guarantee and efficient algorithm, Journal of Machine Learning Research, 22(9):1?32, 2021.
%% [2] Y.C. Yuan, D.F. Sun, and K.C. Toh, An efficient semismooth Newton based algorithm for convex clustering, International Conference on Machine Learning (ICML) 2018.
%%********************************************************************
function [weightVec,NodeArcMatrix,weightMatrix] = compute_weight(X,k,phi,gamma, option)
    if nargin < 5
        option = 0;
    end
    if (gamma <= 0)
       error('Regularized Parameter gamma must be positive.\n');
    end 
    tstart = clock;
    [k_nearest,dist] = knnsearch(X',X','K',k+1);
    len = size(X,2);
    %% construct weight
    weightMatrix = sparse(len,len);
    for i = 1:len
        for j=1:k+1
            weightMatrix(i,k_nearest(i,j)) = exp(-phi*dist(i,j)^2);
            weightMatrix(k_nearest(i,j),i) = exp(-phi*dist(i,j)^2);
        end
        weightMatrix(i,i) = 0;
    end
    %% Temp 
    if option == 1
        for k = 1:1:5
            for i = (k-1)*100+1:1:k*100
                for j = i+1:1:k*100
                    weightMatrix(i,j) = exp(-phi*norm(X(:,i) - X(:,j))^2);
                    weightMatrix(j,i) = weightMatrix(i,j);
                end
            end
        end
    end
    weightMatrix = gamma*weightMatrix;    
    %%construct W, Wbar
    [idx_r,idx_c,val] = find(triu(full(weightMatrix)));
    %[idx_r,idx_c,val] = find(full(weightMatrix));
    num_weight = length(val);
    W = sparse(len,num_weight);
    Wbar = sparse(len,num_weight);
    for i = 1:1:num_weight
        W(idx_r(i),i) = 1;
        Wbar(idx_c(i),i) = 1;
    end
    weightVec = val';
    NodeArcMatrix = W-Wbar;
    fprintf('\ntime taken to generate weight matrix = %3.2f\n',etime(clock,tstart));
%%********************************************************************    
    
    
    
    
    