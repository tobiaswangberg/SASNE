function [A,D,D_knn,ind] = data_to_graph(data,k,knn)
% takes data matrix as input and returns adjacency matrix 
% (0 not connected, 1 connected)
% for connecting the k nearest neighbors.
% if knn is true constructs knn graph, otherwise mutual knn.
    D = squareform(pdist(data));
    D_knn = D;
    [~,ind] = sort(D);
    for i=1:size(D, 1)
        D_knn(i, ind((2 + k):end, i)) = 0;
    end
    if knn
        A = D_knn + D_knn';
        A(A~=0) = 1;
    end
    if ~knn
        A = D .* D';
        A(A~=0) = 1;
    end
    A = sparse(A);
end