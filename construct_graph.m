function [W,A] = construct_graph(data,linear_search,knn,min_k)
    smallest_conn = @smallest_conn;
    D = pdist(data);
    Dsq = squareform(D);
    n = size(Dsq,1);
    [~,ind] = sort(Dsq);
    
    A = smallest_conn(Dsq,linear_search,knn,min_k,'1');
    W = A;
    D_5NN = zeros(n,5);
    for i = 1:n
        D_5NN(i,:) = Dsq(i,ind(2:6,i));
    end

    q = quantile(D_5NN(:),0.9);
    %q = quantile(D(:),0.9);
    Dsq = 6.31*Dsq/q; % normalise to avoid flat region of t-dist
    W(A~=0) = 1./(1 + Dsq(A~= 0).^2);
    W = sparse(W);
 
end

