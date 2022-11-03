%% Compute projection onto the l2 norm ball.
%% Copyright (c) 2021 by
%% Defeng Sun, Kim-Chuan Toh and Yancheng Yuan 
%% Citation:
%% [1] D.F. Sun, K.C. Toh, and Y.C. Yuan, Convex clustering: model, theoretical guarantee and efficient algorithm, Journal of Machine Learning Research, 22(9):1?32, 2021.
%% [2] Y.C. Yuan, D.F. Sun, and K.C. Toh, An efficient semismooth Newton based algorithm for convex clustering, International Conference on Machine Learning (ICML) 2018.
%%*************************************************************************
function output = proj_l2(proj_input, weightVec)
    [d,n] = size(proj_input);
    if n ~= length(weightVec)
       error('Dimensions not agree.\n');
    end
    output = proj_input;
    norm_input = sqrt(sum(proj_input.^2));
    idx = norm_input > weightVec;
     if ~isempty(idx)
        output(:, idx) = bsxfun(@times,output(:,idx),weightVec(idx)./norm_input(idx));   
     end
end
    