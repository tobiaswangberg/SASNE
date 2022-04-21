function pointwise_ranks = dist_to_rank(D)
    D = squareform(D);
    n = size(D,1);
    [~,ind] = sort(D);
    pointwise_ranks = repmat([1:n-1]',1,n);
    for i = 1:n
        r = 1:n;
        r(ind(:,i)) = r;
        r = r(setdiff(1:n,i)) - 1;
        pointwise_ranks(:,i) = r;
    end
end

