function [Y,Z,W] = SASNE(data)
% SASNE - shape aware stochastic neighbour embedding
% Takes as input a HD data matrix (rows are observations, columns features)
% Returns 2D embedding Y, that original graph distance coordinates Z and
% the graph weight matrix W.
    n = size(data,1);
    disp('constructing graph...')
    tic;
    [W,~] = construct_graph(data,true,true,5);
    toc
    t_graph_constr = toc;
    disp('computing graph distance...')
    tic;

    % You can choose graph distance here
    [Z,lambda] = get_symbiharmonic_coords(W,true);
    %[Z,lambda] = get_biharmonic_coords(W,true);
    %[Z,lambda] = get_CTD_coords(W,true);

    toc
    t_graph_dist = toc;
    init_Y = 1e-4.*Z(:,1:2)*sqrt(lambda(2));
    perplexity = floor(0.9*n);
    disp('computing tsne embedding...')
    tic;
    Y = tsne(Z,'InitialY',init_Y,'Exaggeration',12,'LearnRate',...
        n/12,'Perplexity',perplexity,'Verbose',0,'Options',...
        statset('TolFun',1e-100),'Algorithm','exact');
    toc
    t_tsne = toc;
    
    t_total = t_graph_constr + t_graph_dist + t_tsne;
    
    disp(['Total running time ' num2str(t_total) ' seconds']);
end


%%%%%%%%%%%%% HELPER FUNCTION FOR GRAPH DISTANCE %%%%%%%%%%%%%
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
       [V, lambda] = eig(full(Lsym));
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

function [Z,lambdad] = get_CTD_coords(W,exact)
% Computes the commute time distance coordinates
% Takes as input a symmetric weight matrix W where entry (i,j) of W indicates
% the similarity between node i and j along with the boolean argument
% exact. If exact is set to false, then the distance is computed using only
% the 5 % leading eigenvectors to save computational time.
% Returns coordinates Z, the squared Euclidean distances between the column
% vectors of Z is the commute time distance.
    tic
    Lsym = compute_Lsym(W);
    if ~exact
        N = floor(0.05*n);
        [V, lambda] = eigs(sparse(Lsym),N,1e-10);
   end

   if exact
       [V, lambda] = eig(full(Lsym));
   end
 
    [V, lambda] = eig(Lsym);

    % MATLAB function sometimes returns small imaginary part
    V = real(V);
    lambda = real(lambda);
    lambdad = diag(lambda);
    [lambdad,idx] = sort(lambdad,'ascend');
    V = V(:,idx);

    lambdainv = diag(1./sqrt(lambdad(2:end)));
    
    
    deg = sum(W);
    vol = sum(deg);
    sqrt_dinv = 1./sqrt(deg);
    Z = sqrt(vol) * diag(sqrt_dinv)*V(:,2:end)*lambdainv;
    toc

end

function [Z,lambdad] = get_biharmonic_coords(W,exact)
% Computes the biharmonic coordinates
% Takes as input a symmetric weight matrix W where entry (i,j) of W indicates
% the similarity between node i and j along with the boolean argument
% exact. If exact is set to false, then the distance is computed using only
% the 5 % leading eigenvectors to save computational time.
% Returns coordinates Z, the squared Euclidean distances between the column
% vectors of Z is the biharmonic distance distance.
    tic
    d = sum(W,2);
    sqrt_dinv = 1./sqrt(d);


    L = diag(d) - W;
    
    if ~exact
        N = floor(0.05*n);
        [V, lambda] = eigs(sparse(L),N,1e-10);
    end

    if exact
       [V, lambda] = eig(full(L));
    end

    lambdad = diag(lambda);
    [lambdad,idx] = sort(lambdad,'ascend');
    V = V(:,idx);
    lambdainv = diag(1./(lambdad(2:end)));
    disp('Computing biharmonic coordinates...')
    vol = sum(sum(W));

    Z = V(:,2:end)*lambdainv;
    toc

end

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

%%%%%%%%%%%%% HELPER FUNCTION FOR GRAPH CONSTRUCTION %%%%%%%%%%%%%
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
        D_5NN(i,:) = Dsq(i,ind(2:6,i));
    end

    q = quantile(D_5NN(:),0.9);
    Dsq = 6.31*Dsq/q; % normalise to avoid flat region of t-dist
    
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

function A = smallest_conn(D,linear_search,knn,min_k,layer)
    if linear_search
        smallest_k = @smallest_k_linear;
    end
    if ~linear_search
        smallest_k = @smallest_k_binary;
    end
   
    [~,ind] = sort(D);

    k = smallest_k(D,ind,min_k,knn);
    A = DtoA(D,k,ind,knn);
    %disp(['Smallest k is ',int2str(k),' in layer ',layer])
    if k > min_k
        k_disconnect = k - 1;
        A_disconnect = DtoA(D,k_disconnect,ind,knn);
        [comps,N] = find_comps(A_disconnect);
        disp(['we have ',int2str(N),' connected components for k-1 = ',int2str(k-1),' in layer ',layer])
        for i = 1:N
            disp(['Looking into component ',int2str(i)])
            comp_ind = find(comps == i);
            sub_D = D(comp_ind,comp_ind);
            new_layer = [layer,'.',int2str(i)];
            A_sub = smallest_conn(sub_D,linear_search,knn,min_k,new_layer);
            A(comp_ind,comp_ind) = A_sub;
        end
    end
    
end

function A = DtoA(D,k,ind,knn)
    for i=1:size(D, 1)
        D(i, ind((2 + k):end, i)) = 0;
    end
    if knn
        A = D + D';
        A(A~=0) = 1;
    end
    if ~knn
        A = D .* D';
    end
    A = sparse(A);
end

function k = smallest_k_binary(D,ind,min_k,knn)
    n = size(D,1);    
    k_lower = min_k;
    k_upper = ceil(n/2);
    k = floor((k_lower + k_upper)/2);
    
    connected = is_conn(D,ind,k,knn);
    while (k_upper - k_lower) > 1
        connected = is_conn(D,ind,k,knn);    
        if connected
            %disp('connected')
            k_upper = k;
            k = floor((k_lower + k_upper)/2);
        end
        if ~connected
            %disp('not connected')
            k_lower = k;
            k = ceil((k_lower + k_upper)/2);
        end
    end  
    k = max(k_upper,min_k);
end

function k = smallest_k_linear(D,ind,min_k,knn)   
    k = min_k;
    
    connected = is_conn(D,ind,k,knn);
    while ~connected
        k = k + 1;
        connected = is_conn(D,ind,k,knn);
    end

end


function connected = is_conn(D,ind,k,knn)
% Returns true if graph is connected
    A = DtoA(D,k,ind,knn);
    connected = is_connected(A);
end

function connected = is_connected(A)
    n = size(A,1);
    % initiliase the recursive dfs
    marked = 1; 
    unmarked = 2:n;
    v_curr = 1;

    [comp,marked,unmarked] = dfs(v_curr,marked,unmarked,A);
    connected = isempty(unmarked);
end

function [comps,count] = find_comps(A)
% Takes as input an adjacency matrix A
% returns the disconneted components of the graph comps and the number
% of disconnected components.
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

function [comp,marked,unmarked] = dfs(v,marked,unmarked,A)
    NNs = A(v,:);
    NNs = find(NNs > 0);
    %disp(['knns ',int2str(length(knns)),newline])
    NNs = intersect(NNs,unmarked);
    comp = NNs;
    marked = union(marked,comp);
    unmarked = setdiff(unmarked,comp);
    if ~isempty(comp)
        n = length(comp);
        for i = 1:n
            [new_comp,marked,unmarked] = dfs(NNs(i),marked,unmarked,A);
            comp = union(comp,new_comp);
        end
    end
end

