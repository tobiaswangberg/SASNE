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