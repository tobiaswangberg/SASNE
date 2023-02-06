function A = construct_graph_fast(data,k)
% construct graph by connecting the k nearest neighbors
% if graph is not connected, locate the disconnected components
% and then connect the k nearest nodes between the component and 
% the rest of the graph
    [A,D] = data_to_graph(data,k,true);
    connected = is_connected(A);
    if ~connected
        disp(['graph is not connected for  k = ',num2str(k),newline,...
            'connecting disconnected components'])
        [comps,N] = find_comps(A,k);
        for i = 1:N-1
            for j = (i+1):N
                comp_ind = find(comps == i);
                other_comp_ind = find(comps == j);
                D_sub = D(comp_ind,other_comp_ind);
                [nrow,ncol] = size(D_sub);
                [row,col,vals] = mink_mat(D_sub,k);
               
                for l = 1:k
                    A(comp_ind(row(l)),other_comp_ind(col(l))) = 1;
                    A(other_comp_ind(col(l)),comp_ind(row(l))) = 1;
                end
%                 A(comp_ind(row),other_comp_ind(col))...
%                     = A(comp_ind(row),other_comp_ind(col)) + eye(k);
%                 A(other_comp_ind(col),comp_ind(row))...
%                     = A(other_comp_ind(col),comp_ind(row)) + eye(k);
%                 A = A + A';
                %A(A~=0) = 1;
            end
        end
    end
end


function connected = is_conn(D,ind,k,knn)
% Returns true if graph is connected
    A = DtoA(D,k,ind,knn);
    connected = is_connected(A);
end

function [connected,comp,marked,unmarked] = is_connected(A)
    n = size(A,1);
    % initiliase the recursive dfs
    marked = 1; 
    unmarked = 2:n;
    v_curr = 1;
    [comp,marked,unmarked] = dfs(v_curr,marked,unmarked,A);
    connected = isempty(unmarked);
end

function [comps,count] = find_comps(A,k)
    n = size(A,1);
    unmarked = 1:n;
    count = 0;
    comps = ones(1,n);
    while ~isempty(unmarked)
        marked = [];
        v_curr = unmarked(1);
        marked = union(v_curr,marked);
        unmarked = setdiff(unmarked,v_curr);
        [comp,marked,unmarked] = dfs(v_curr,marked,unmarked,A);
        count = count + 1;
        comps(marked) = count*ones(length(marked),1);
    end

end


function [comp,marked,unmarked] = dfs(v,marked,unmarked,A)
    NNs = A(v,:);
    NNs = find(NNs > 0);
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
