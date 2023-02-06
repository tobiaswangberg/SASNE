function k = smallest_k_linear(D,ind,min_k,knn)   
    k = min_k;
    
    connected = is_conn(D,ind,k,knn);
    while ~connected
        k = k + 1;
        connected = is_conn(D,ind,k,knn);
    end

end