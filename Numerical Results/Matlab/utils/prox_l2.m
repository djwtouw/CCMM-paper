%%********************************************************************
%% Compute proximal mapping for l2 norm.
%% Copyright (c) 2021 by
%% Defeng Sun, Kim-Chuan Toh and Yancheng Yuan 
%% Citation:
%% [1] D.F. Sun, K.C. Toh, and Y.C. Yuan, Convex clustering: model, theoretical guarantee and efficient algorithm, Journal of Machine Learning Research, 22(9):1?32, 2021.
%% [2] Y.C. Yuan, D.F. Sun, and K.C. Toh, An efficient semismooth Newton based algorithm for convex clustering, International Conference on Machine Learning (ICML) 2018.
%%*************************************************************************
function [x, rr, norm_y_col] = prox_l2(y,w)
    [d,n] = size(y);
    if n ~= length(w)
       error('Dimensions not agree.\n');
    end
    norm_y_col = sqrt(sum(y.*y));
    alpha_vec = w./(norm_y_col + 1e-15);
    rr = alpha_vec < 1;
    idx = find(rr); 
    x = sparse(d,n);
    if ~isempty(idx)
        x(:,idx) = bsxfun(@times,y(:,idx),1-alpha_vec(idx)); 
    end
%%********************************************************************    
    