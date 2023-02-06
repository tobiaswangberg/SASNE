function connected = is_conn(D,ind,k,knn)
% Returns true if graph is connected
    A = DtoA(D,k,ind,knn);
    connected = is_connected(A);
end