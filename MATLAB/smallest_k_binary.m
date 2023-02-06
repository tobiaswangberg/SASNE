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