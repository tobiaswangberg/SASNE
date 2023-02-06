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