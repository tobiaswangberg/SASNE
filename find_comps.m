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
        [marked,unmarked] = dfs(v_curr,marked,unmarked,A);
        count = count + 1;
        comps(marked) = count*ones(length(marked),1);
    end

end



function [marked,unmarked] = dfs(v,marked,unmarked,A)
    NNs = A(v,:);
    NNs = find(NNs > 0);
    %disp(['knns ',int2str(length(knns)),newline])
    NNs = intersect(NNs,unmarked);

    marked = union(marked,NNs);
    unmarked = setdiff(unmarked,marked);
    if ~isempty(NNs)
        n = length(NNs);
        for i = 1:n
            [new_marked,unmarked] = dfs(NNs(i),marked,unmarked,A);
            marked = union(marked,new_marked);
        end
    end
end