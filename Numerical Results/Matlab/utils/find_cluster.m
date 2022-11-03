%% Visualization of the clustering results
%%-----------------------------------------------------------------
%% X   = each column is the computed centroid associated 
%%       with the corresponding data point.
%% tol = tolerance to declare that two centroids are identical.
%%-----------------------------------------------------------------
    function [cluster_id, num_cluster] = find_cluster(X,tol)
    
    if isempty(X)
       error('Input Data must be nonempty.\n');
    end
    [m,n] = size(X);
    index_temp = [1:n];
    size_index = n;
    cluster_id = ones(n,1);
    reference_id = 0;
    while size_index > 0
        reference_id = reference_id + 1;
        reference_point = X(:,index_temp(1));
        cluster_id(index_temp(1)) = reference_id;
        index_temp = index_temp(2:end); %setdiff(index_temp,index_temp(1));
        if (true) 
           Xdiff = X(:,index_temp)-reference_point*ones(1,length(index_temp));
           normXdiff = sqrt(sum(Xdiff.*Xdiff));
           idx = find(normXdiff < tol*m);
           index_t = index_temp(idx); 
           cluster_id(index_t) = reference_id;
        end
        index_temp = setdiff(index_temp,index_t);
        size_index = length(index_temp);
    end
    num_cluster = reference_id;
end
    