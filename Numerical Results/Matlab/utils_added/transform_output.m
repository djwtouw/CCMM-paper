%**************************************************************************
% Transforms output of the SSNAL algorithm into the clusterpath format used
% by the CCMM algorithm
%**************************************************************************

function path = transform_output(ssnal_object)
    % Get dimensions of the data
    n = ssnal_object{1, 1}.n;
    d = ssnal_object{1, 1}.d;

    % Get the number of gammas used for the clusterpath
    n_gammas = size(ssnal_object);
    n_gammas = n_gammas(2);

    % Initialize the result matrix
    path = zeros(n * n_gammas, d + 2);

    % Fill in the results iteratively for each gamma
    for i = 1:n_gammas
        path((n * (i - 1) + 1):(n * i), 1) = 1:n;
        path((n * (i - 1) + 1):(n * i), 2) = ssnal_object{1, i}.gamma;
        path((n * (i - 1) + 1):(n * i), 3:(d + 2)) = ssnal_object{1, i}.solution.X';
    end
end
