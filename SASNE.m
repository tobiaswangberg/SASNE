function Y = SASNE(data)
% SASNE - shape aware stochastic neighbour embedding
% Takes as input a HD data matrix (rows are observations, columns features)
% Returns 2D embedding Y
    n = size(data,1);
    construct_graph = @construct_graph;
    disp('constructing graph...')
    tic;
    [W,~] = construct_graph(data,true,true,5);
    t_graph = toc;
    tic;
    [Z,lambda] = get_symbiharmonic_coords(W,true);
    clear W
    t_diag = toc;
    init_Y = 1e-4.*Z(:,1:2)*sqrt(lambda(2));
    perplexity = floor(0.9*n);
    disp('Computing tsne embedding...')
    tic;
    Y = tsne(Z,'InitialY',init_Y,'Exaggeration',12,'LearnRate',...
        n/12,'Perplexity',perplexity,'Verbose',0,'Options',...
        statset('TolFun',1e-100),'Algorithm','exact');
    t_tsne = toc;
    
    t_total = t_graph + t_diag + t_tsne;
    
    disp(['Total time ' num2str(t_total) ' seconds']);
end

