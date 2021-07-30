function [A] = smallest_connK_discrete(D,upper_k,knn,k_start)

if knn
    sym_func = @sym;
    disp('construncting knn adjacency matrix...')
else
    sym_func = @mSym;
    disp('constructing mutual knn adjacency matrix...')
end

k=k_start;

[tmp,ind] = sort(D);
n = size(D,1);
tic
disp('Finding smallest connected k...')
[connected,A] = compute_A(D,ind,k,knn);

while ~connected
    %disp(['Graph still not connected for k = ',num2str(k),'... increasing k by 1'])
    k = k+1;
    [connected,A] = compute_A(D,ind,k,knn);
    if k == upper_k
        disp('Upper bound on k reached, connecting disconnected components...')
        A = connect_comps(A,D);
        break
    end
end
disp(['Graph is connected for k is equal to ',num2str(k)])

toc

end

function [connected,G] = compute_A(D,ind,k,knn)
if knn
    sym_func = @sym;
else
    sym_func = @mSym;
end
n = size(D,1);
for i=1:size(D, 1)
    D(i, ind((2 + k):end, i)) = 0;
end

D(D ~= 0) = 1;
G = sym_func(D);
d = sum(G,1);
if sum(d==0) > 0 % prevent zero degree outliers by connecting to nearest neighbor
    
    ind0 = find(d==0);
    nr_0 = max(size(ind0));
    disp(['Placing edge to ',num2str(nr_0),' isolated vertices...'])
    
    for i = 1:nr_0
        G(ind0(i),ind(2:7,ind0(i))) = 1;
        G(ind(2:7,ind0(i)),ind0(i)) = 1;
    end
    d = sum(G,2);
end
I = eye(size(G,1));

dinv = 1./sqrt(d);

W = G;
for i = 1:n
    for j = 1:n
        W(i,j) = W(i,j)*(dinv(i)*dinv(j));
    end
end

Lsym = I - W;

lambdas = eigs(sparse(Lsym),10,1e-15);
nr_disconnected = sum(lambdas(2:end) < 1e-5);
connected = lambdas(2) > 1e-5;
if ~connected
    disp([newline,'Graph still not connected for k = ',num2str(k),'...,',newline,...
        'nr of disconnected components are ',(num2str(nr_disconnected+1))]);
end
end

function Gsym = mSym(G)
Gsym = G.*G';
end
function Gsym = sym(G)
Gsym = G + G';
Gsym(Gsym~=0) = 1;
end

function connect_A = connect_comps(A,D)

[foo, p, bar, r] = dmperm(A);
sizes = diff(r)
k = length(sizes)
% Now compute the array blocks
n = size(A,1);
blocks = zeros(1, n);
blocks(r(1:k)) = ones(1, k);
blocks = cumsum(blocks);

% Permute blocks so it maps vertices of A to components
blocks(p) = blocks;
for i = 1:k
    ind_i = find(blocks == i);
    not_i = setdiff(1:k,i);
    for j = 1:k-1
        curr_j = not_i(j);
        ind_j = find(blocks == curr_j);
        D_curr = D(ind_i,ind_j);
        
        nr_connects = min(size(D_curr,1),5);
        k_smallest = mink(D_curr(D_curr>0),nr_connects);
        [x,y] = find(D_curr >= k_smallest(1) & D_curr <= k_smallest(nr_connects));
        A(ind_i(x),ind_j(y)) = 1;
        A(ind_j(y),ind_i(x))= 1;
    end
end
d = sum(A,2);
connect_A = A;

end