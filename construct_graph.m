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

function [comps,count] = find_comps(A,k)
    n = size(A,1);

    
    unmarked = 1:n;
    count = 0;
    comps = ones(1,n);
    while ~isempty(unmarked)
        marked = [];
        v_curr = unmarked(1);
        marked = union(v_curr,marked);
        unmarked = setdiff(unmarked,v_curr);
        [marked,unmarked] = dfs(v_curr,marked,unmarked,A);
        count = count + 1;
        comps(marked) = count*ones(length(marked),1);
    end

end



function [marked,unmarked] = dfs(v,marked,unmarked,A)
    NNs = A(v,:);
    NNs = find(NNs > 0);
    %disp(['knns ',int2str(length(knns)),newline])
    NNs = intersect(NNs,unmarked);

    marked = union(marked,NNs);
    unmarked = setdiff(unmarked,marked);
    if ~isempty(NNs)
        n = length(NNs);
        for i = 1:n
            [new_marked,unmarked] = dfs(NNs(i),marked,unmarked,A);
            marked = union(marked,new_marked);
        end
    end
end
