function Lsym = compute_Lsym(W)
%%% Function takes as input symmetric weight matrix W and returns symmetric
%%% symmetric Laplacian Lsym 
    d = sum(W,2);

    sqrt_dinv = 1./sqrt(d);
    n = size(W,1);
    % compute normalised weight matrix make use of sparsity to save some
    % time
    Wnorm = zeros(n,n); 
    nnzs = nnz(W);
    [r,c] = find(W);
 
    for i = 1:nnzs
            Wnorm(r(i),c(i)) = W(r(i),c(i)) * sqrt_dinv(r(i))*sqrt_dinv(c(i));
    end

    I = eye(size(Wnorm,1));
    Lsym = I - Wnorm;
end