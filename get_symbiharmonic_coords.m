function [Z,lambdad] = get_symbiharmonic_coords(W,exact)
    disp('computing symmetric biharmonic coordinates')
    tic
    if ~issparse(W)
        W = sparse(W);
    end
  

    d = sum(W,2);

    sqrt_dinv = 1./sqrt(d);
    n = size(W,1);
    Wnorm = zeros(n,n); % compute normalised weight matrix
    nnzs = nnz(W);
    [r,c] = find(W);
 
    for i = 1:nnzs
            Wnorm(r(i),c(i)) = W(r(i),c(i)) * sqrt_dinv(r(i))*sqrt_dinv(c(i));
    end

    I = eye(size(Wnorm,1));
    Lsym = I - Wnorm;
    disp('computing eigenspectrum')
    
    if ~exact
        sum_lambda = trace(Lsym);
        N = floor(0.05*n);
        [V, lambda] = eigs(sparse(Lsym),N,1e-10);
        disp([num2str(sum(lambda(:))/sum_lambda),' % of variance explained'])
    end

   if exact
       [V, lambda] = eig(full(Lsym));
   end
       
    lambdad = diag(lambda);
    [lambdad,idx] = sort(lambdad,'ascend');
    V = V(:,idx);

    lambdainv = sparse(diag(1./(lambdad(2:end))));
    
    vol = sum(sum(W));


    Z = sqrt(vol) * sparse(diag(sqrt_dinv))*V(:,2:end)*lambdainv;
    toc

end
