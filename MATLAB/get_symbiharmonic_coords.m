function [Z,lambdad] = get_symbiharmonic_coords(W,exact)
% Computes the symmetric biharmonic graph distance coordinates
% Takes as input a symmetric weight matrix W where entry (i,j) of W indicates
% the similarity between node i and j along with the boolean argument
% exact. If exact is set to false, then the distance is computed using only
% the 5 % leading eigenvectors to save computational time.
% Returns coordinates Z, the squared Euclidean distances between the column
% vectors of Z is the biharmonic distance, computed using the symmetric 
% Laplacian.
    if ~issparse(W)
        W = sparse(W); 
    end
    
    Lsym = compute_Lsym(W);
    if Lsym ~= Lsym'
        disp('Warning symmetric Laplacian is not symmetric!')
    end
   
   
    if ~exact
        sum_lambda = trace(Lsym);
        N = floor(0.05*n);
        [V, lambda] = eigs(sparse(Lsym),N,1e-10);
        disp([num2str(sum(lambda(:))/sum_lambda),' % of variance explained'])
    end

   if exact
       [V, lambda] = eig(full(Lsym))
   end
   % MATLAB function sometimes returns small imaginary part for numerical
   % reasons
    V = real(V);
    lambda = real(lambda);  
    lambdad = diag(lambda);
    
    [lambdad,idx] = sort(lambdad,'ascend');
    V = V(:,idx);
    
    lambdainv = sparse(diag(1./(lambdad(2:end))));
    
    deg = sum(W);
    vol = sum(deg);
    sqrt_dinv = 1./sqrt(deg);

    
    Z = sqrt(vol) * sparse(diag(sqrt_dinv))*V(:,2:end)*lambdainv;

end