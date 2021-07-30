function [Z,lambdad,d,vol,W] = compute_Z(data,upper_k,knn,start_k)
if ~exist('knn', 'var')
    knn = true;
end
if ~exist('start_k', 'var')
    start_k = 2;
end
if ~exist('upper_k', 'var')
    upper_k = size(data,1);
end
disp('Computing Euclidean distance matrix...')
tic
D = squareform(pdist(data));
toc
A = smallest_connK_discrete(D,upper_k,knn,start_k); % return adjacency matrix
tic
disp('Computing eigenspectrum of symmetric Laplacian')
W = D;
W(A~=0) = D(A~=0).^(-2); % compute weighted adjacency matrix
W(A==0) = 0;
d = sum(W,2);
sqrt_dinv = 1./sqrt(d);

Wnorm = W; % compute normalised weight matrix
for i = 1:size(Wnorm,1)
    for j=1:size(Wnorm,2)
        Wnorm(i,j) = Wnorm(i,j) .* (sqrt_dinv(i)*sqrt_dinv(j));
    end
end



I = eye(size(Wnorm,1));
Lsym = I - Wnorm;
[V, lambda] = eig(full(Lsym));
lambdad = diag(lambda);



if lambdad(2) < 1e-5
    disp('Warning: Graph does not appear to be fully connected.')
end
toc

tic
disp('Computing CTD coordinates...')
vol = sum(sum(W));

Z = sqrt(vol) * diag(1./sqrt(d))*V(:,2:end)*diag(1./lambdad(2:end));
toc

end



