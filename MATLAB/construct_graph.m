function [W,A] = construct_graph(data,linear_search,knn,min_k)
 % Function to create a graph representation of the data.
 % The function takes as input a data matrix (rows observations, colums
 % features), 
 % boolean variable linear search. If this argument is set to true,
 % the search for the smallest k such that the graph is connected is done 
 % linearly by increasing k starting from a low value. Otherwise this is
 % done by a binary search.
 % boolean argument knn. If set to true then the graph is constructed by
 % connecting the k nearest neighbors, otherwise a mutual kNN graph is
 % constructured.
 % integer argument min_k indicates the minimum numbers k nearest neighbors
 % of the graph, so that the graph construction guarantees at least min_k 
 % neighbors of each node are connected in the graph. 
 % Returns the weight matrix W of the graph and binary adjacency matrix A.

    % need to handle case when discrete data is not unique to avoid
    % degenerate case when distance would be 0 between points.
    [data_unique,ia,ic] = unique(data, 'rows', 'stable');
    Cnts = accumarray(ic, 1);

    D = pdist(data_unique);
    Dsq = squareform(D);
    n = size(Dsq,1);
    [~,ind] = sort(Dsq);
    
    A = smallest_conn(Dsq,linear_search,knn,min_k,'1');
    W = A;
    D_5NN = zeros(n,min_k);
    for i = 1:n
        D_5NN(i,:) = Dsq(i,ind(2:(min_k+1),i));
    end

    q = quantile(D_5NN(:),0.9);
    %Dsq = 6.31*Dsq/q; % normalise to avoid flat region of t-dist
    
    nnzs = nnz(A);
    [r,c] = find(A);
    W = A;
 
    for i = 1:nnzs
            W(r(i),c(i)) = 1/(1 + Dsq(r(i),c(i))^2)*Cnts(r(i))*Cnts(c(i));
            W(c(i),r(i)) = W(r(i),c(i));
    end
    %W(A~=0) = 1./(1 + Dsq(A~= 0).^2);
    W = sparse(W);
 
end